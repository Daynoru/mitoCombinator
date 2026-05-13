use anyhow::Result;
use fxhash::FxHashMap;
use needletail::{parse_fastx_file, Sequence};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Write, stdout};
use std::path::Path;
use std::time::Instant;

mod restr;

const K: usize = 12;
const MIN_REPEAT: usize = 50;
const MAX_REPEAT: usize = 1500;
const MIN_DISTANCE: usize = 5000;
const FLANK: usize = 20000;

fn progress(step: usize, total: usize, label: &str) {
    let w = 40;
    let pct = step * 100 / total.max(1);
    let filled = step * w / total.max(1);
    let bar: String = (0..w).map(|i| if i < filled { '█' } else { '░' }).collect();
    print!("\r\x1b[K  {} [{bar}] {:>3}% ({}/{})", label, pct, step, total);
    stdout().flush().ok();
}

fn main() -> Result<()> {
    let start = Instant::now();

    let fasta_path = find_fasta("data")?;
    let mut reader = parse_fastx_file(&fasta_path)?;
    let record = reader.next().ok_or_else(|| anyhow::anyhow!("Empty FASTA"))??;
    let seq = record.sequence();
    let seq_upper: Vec<u8> = seq
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        .collect();
    let seq_len = seq_upper.len();

    // --- Stage 1: k-mer index ---
    eprint!("[1/5] Building k-mer index... ");
    let mut kmer_map: FxHashMap<u64, Vec<usize>> = FxHashMap::default();
    for i in 0..seq_len.saturating_sub(K) {
        kmer_map.entry(hash_kmer(&seq_upper[i..i + K])).or_default().push(i);
    }
    eprintln!("done");

    // --- Stage 2: repeat detection ---
    eprintln!("[2/5] Detecting repeats...");
    let mut repeats: Vec<(usize, usize, usize, usize)> = Vec::new();
    let total_hashes = kmer_map.len();
    let mut processed = 0;
    for positions in kmer_map.values() {
        processed += 1;
        if processed % 50000 == 0 {
            progress(processed, total_hashes, "  k-mer groups");
        }
        if positions.len() < 2 { continue; }
        let mut by_kmer: FxHashMap<Vec<u8>, Vec<usize>> = FxHashMap::default();
        for &pos in positions {
            by_kmer.entry(seq_upper[pos..pos + K].to_vec()).or_default().push(pos);
        }
        for grouped in by_kmer.values() {
            if grouped.len() < 2 { continue; }
            for i in 0..grouped.len() {
                for j in (i + 1)..grouped.len() {
                    let p1 = grouped[i];
                    let p2 = grouped[j];
                    let mut left = 0usize;
                    while left < p1.min(p2) && seq_upper[p1 - left - 1] == seq_upper[p2 - left - 1] { left += 1; }
                    let mut right = K;
                    while p1 + right < seq_len && p2 + right < seq_len && seq_upper[p1 + right] == seq_upper[p2 + right] { right += 1; }
                    let len = left + right;
                    if len >= MIN_REPEAT && len <= MAX_REPEAT {
                        let s1 = p1 - left;
                        let s2 = p2 - left;
                        repeats.push((s1, s1 + len, s2, s2 + len));
                    }
                }
            }
        }
    }
    progress(total_hashes, total_hashes, "  k-mer groups");
    eprintln!();
    repeats.sort_by_key(|r| (r.0, r.2));
    let mut seen = HashSet::new();
    repeats.retain(|r| seen.insert((r.0, r.1, r.2, r.3)));
    let filtered: Vec<_> = repeats.into_iter()
        .filter(|&(s1, e1, s2, e2)| distance(s1, e1, s2, e2, seq_len) >= MIN_DISTANCE)
        .collect();
    eprintln!("  {} direct repeats found", filtered.len());

    // --- Export repeats ---
    let mut rep_file = File::create("repeats.csv")?;
    writeln!(rep_file, "id,start1,end1,start2,end2,length,distance")?;
    for (i, &(s1, e1, s2, e2)) in filtered.iter().enumerate() {
        writeln!(rep_file, "{},{},{},{},{},{},{}", i + 1, s1, e1, s2, e2, e1 - s1, distance(s1, e1, s2, e2, seq_len))?;
    }

    // --- Stage 3: restriction map ---
    eprint!("[3/5] Building restriction map... ");
    let enzymes = restr::parse_enzymes("data/enzymes.csv")?;
    let mut rmap = restr::build_map(seq, &enzymes);
    rmap.retain(|_, s| s.len() <= 500);
    let mut nmap: HashMap<String, Vec<usize>> = HashMap::new();
    for (name, sites) in &rmap {
        nmap.entry(normalize_name(name)).or_default().extend(sites);
    }
    eprintln!("{} enzymes", nmap.len());

    // --- Stage 4: diagnostic fragments ---
    eprintln!("[4/5] Searching diagnostic fragments...");
    let mut all_sol: Vec<(usize, String, String, char, usize)> = Vec::new();
    let total_repeats = filtered.len();
    for (idx, &(s1, e1, s2, e2)) in filtered.iter().enumerate() {
        if idx % 50 == 0 {
            progress(idx, total_repeats, "  repeats");
        }
        let r = Repeat { start1: s1, end1: e1, start2: s2, end2: e2, length: e1 - s1 };
        for sol in find_diag(&r, &nmap, seq_len) {
            all_sol.push((idx, sol.enzyme1, sol.enzyme2, sol.ring, sol.fragment_length));
        }
    }
    progress(total_repeats, total_repeats, "  repeats");
    eprintln!();
    eprintln!("  {} diagnostic variants", all_sol.len());

    // --- Stage 5: optimization ---
    eprint!("[5/5] Optimization... ");
    let mut cov: HashMap<(String, String), HashSet<usize>> = HashMap::new();
    for (idx, e1, e2, _, _) in &all_sol {
        let key = if e1 <= e2 { (e1.clone(), e2.clone()) } else { (e2.clone(), e1.clone()) };
        cov.entry(key).or_default().insert(*idx);
    }
    let mut ranked: Vec<_> = cov.iter().map(|((e1, e2), r)| (e1.clone(), e2.clone(), r.len())).collect();
    ranked.sort_by_key(|(_, _, c)| -(*c as i64));
    let selected = greedy_select(&cov, filtered.len());
    let mut cum: HashSet<usize> = HashSet::new();
    eprintln!("done\n");

    // --- Output ---
    println!("Top-10 enzyme pairs:");
    for (e1, e2, c) in ranked.iter().take(10) {
        println!("  {} + {} : {} repeats", e1, e2, c);
    }
    println!("Optimal set:");
    for (i, (e1, e2)) in selected.iter().enumerate() {
        if let Some(r) = cov.get(&(e1.clone(), e2.clone())) {
            cum.extend(r);
            println!("  {}. {} + {} : {} covered (cumulative {}/{})", i + 1, e1, e2, r.len(), cum.len(), filtered.len());
        }
    }

    // --- Export diagnostics ---
    let mut f = File::create("diagnostic_results.csv")?;
    writeln!(f, "repeat_idx,enzyme1,enzyme2,ring,fragment_length,in_optimal_set")?;
    let opt: HashSet<(String, String)> = selected.iter().cloned().collect();
    for (idx, e1, e2, ring, len) in all_sol.iter().take(10000) {
        let key = if e1 <= e2 { (e1.clone(), e2.clone()) } else { (e2.clone(), e1.clone()) };
        writeln!(f, "{},{},{},{},{},{}", idx, e1, e2, ring, len, opt.contains(&key))?;
    }
    eprintln!("Done in {:.2} s", start.elapsed().as_secs_f64());
    Ok(())
}

// ===== Data structures =====
struct Repeat { start1: usize, end1: usize, start2: usize, end2: usize, length: usize }
struct DiagSol { enzyme1: String, enzyme2: String, ring: char, fragment_length: usize }

// ===== Distance (circular) =====
fn distance(s1: usize, e1: usize, s2: usize, e2: usize, g: usize) -> usize {
    let d1 = if s2 >= e1 { s2 - e1 } else { g - e1 + s2 };
    let d2 = if s1 >= e2 { s1 - e2 } else { g - e2 + s1 };
    d1.min(d2)
}

// ===== K-mer hash =====
fn hash_kmer(kmer: &[u8]) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    for &b in kmer.iter() { b.to_ascii_uppercase().hash(&mut hasher); }
    hasher.finish()
}

// ===== Enzyme name normalization =====
fn normalize_name(name: &str) -> String {
    name.replace("WarmStart ", "").replace(" §", "").replace("§", "")
        .trim_end_matches("-HFv2").trim_end_matches("-HF").trim_end_matches("-v2")
        .trim().to_string()
}

// ===== Greedy set cover =====
fn greedy_select(cov: &HashMap<(String, String), HashSet<usize>>, total: usize) -> Vec<(String, String)> {
    let mut sel = Vec::new();
    let mut covered: HashSet<usize> = HashSet::new();
    let mut pairs: Vec<_> = cov.iter().collect();
    pairs.sort_by_key(|(_, r)| -(r.len() as i64));
    for _ in 0..pairs.len() {
        if covered.len() >= total { break; }
        let mut best: Option<((String, String), usize)> = None;
        for ((e1, e2), reps) in &pairs {
            let new = reps.difference(&covered).count();
            if new > 0 && best.as_ref().map_or(true, |(_, c)| new > *c) {
                best = Some(((e1.clone(), e2.clone()), new));
            }
        }
        if let Some(((e1, e2), _)) = best {
            sel.push((e1.clone(), e2.clone()));
            if let Some(r) = cov.get(&(e1.clone(), e2.clone())) { covered.extend(r); }
        } else { break; }
    }
    sel
}

// ===== Diagnostic fragment search =====
fn find_diag(r: &Repeat, rmap: &HashMap<String, Vec<usize>>, g: usize) -> Vec<DiagSol> {
    let mut sol = Vec::new();
    let f = FLANK;
    let r1 = r.start1.saturating_sub(f)..=(r.start1 + f).min(g);
    let r2 = r.end1.saturating_sub(f)..=(r.end1 + f).min(g);
    let r3 = r.start2.saturating_sub(f)..=(r.start2 + f).min(g);
    let r4 = r.end2.saturating_sub(f)..=(r.end2 + f).min(g);
    for (e1, s1) in rmap {
        let a1: Vec<usize> = s1.iter().filter(|&&p| r1.contains(&p)).copied().collect();
        let a2: Vec<usize> = s1.iter().filter(|&&p| r2.contains(&p)).copied().collect();
        for (e2, s2) in rmap {
            let b1: Vec<usize> = s2.iter().filter(|&&p| r3.contains(&p)).copied().collect();
            let b2: Vec<usize> = s2.iter().filter(|&&p| r4.contains(&p)).copied().collect();
            for &p1 in &a1 { for &p2 in &a2 { for &p3 in &b1 { for &p4 in &b2 {
                let nr1 = p2.saturating_sub(p1);
                let nr2 = p4.saturating_sub(p3);
                let rr1 = p4.saturating_sub(p1);
                let rr2 = p2.saturating_sub(p3);
                if (nr1 != rr1 || nr2 != rr2) && ((rr1 >= 500 && rr1 <= 10000) || (rr2 >= 500 && rr2 <= 10000)) {
                    sol.push(DiagSol { enzyme1: e1.clone(), enzyme2: e2.clone(), ring: 'A', fragment_length: rr1.max(rr2) });
                }
            }}}}
        }
    }
    sol
}

fn find_fasta(dir: &str) -> Result<std::path::PathBuf> {
    let dir = Path::new(dir);
    if !dir.exists() { anyhow::bail!("Directory '{}' not found", dir.display()); }
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if let Some(ext) = path.extension() {
            if ext == "fasta" || ext == "fa" || ext == "fna" { return Ok(path); }
        }
    }
    anyhow::bail!("No .fasta file found in '{}'", dir.display())
}
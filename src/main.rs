use anyhow::Result;
use fxhash::FxHashMap;
use needletail::{parse_fastx_file, Sequence};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

mod restr;

const K: usize = 15;
const MIN_REPEAT: usize = 50;
const MAX_REPEAT: usize = 1500;
const MIN_DISTANCE: usize = 5000;
const FLANK: usize = 3500;

fn main() -> Result<()> {
    let start = Instant::now();

    let fasta_path = find_fasta("data")?;
    println!("Референс: {}", fasta_path.display());

    let mut reader = parse_fastx_file(&fasta_path)?;
    let record = reader
        .next()
        .ok_or_else(|| anyhow::anyhow!("FASTA файл пуст"))??;
    let seq = record.sequence();
    let seq_len = seq.len();
    println!("Длина генома: {} п.н.", seq_len);

    // --- Этап 1: поиск повторов ---
    println!("Строим хеш-таблицу k-меров (k={})...", K);
    let mut kmer_map: FxHashMap<u64, Vec<usize>> = FxHashMap::default();

    for i in 0..seq_len.saturating_sub(K) {
        let kmer = &seq[i..i + K];
        let hash = hash_kmer(kmer);
        kmer_map.entry(hash).or_default().push(i);
    }
    println!(
        "Уникальных k-меров: {} (из {} всего)",
        kmer_map.len(),
        seq_len - K + 1
    );

    println!("Ищем повторы...");
    let mut repeats: Vec<Repeat> = Vec::new();

    for (_hash, positions) in &kmer_map {
        if positions.len() < 2 {
            continue;
        }
        for i in 0..positions.len() {
            for j in (i + 1)..positions.len() {
                let pos1 = positions[i];
                let pos2 = positions[j];

                let mut left = 0;
                while left < pos1.min(pos2)
                    && seq[pos1 - left - 1].to_ascii_uppercase()
                        == seq[pos2 - left - 1].to_ascii_uppercase()
                {
                    left += 1;
                }

                let mut right = K;
                while pos1 + right < seq_len
                    && pos2 + right < seq_len
                    && seq[pos1 + right].to_ascii_uppercase()
                        == seq[pos2 + right].to_ascii_uppercase()
                {
                    right += 1;
                }

                let total_len = left + right;
                let start1 = pos1 - left;
                let start2 = pos2 - left;

                if total_len >= MIN_REPEAT && total_len <= MAX_REPEAT {
                    repeats.push(Repeat {
                        start1,
                        end1: start1 + total_len,
                        start2,
                        end2: start2 + total_len,
                        length: total_len,
                    });
                }
            }
        }
    }

    println!("Всего пар до фильтрации: {}", repeats.len());
    let filtered = filter_repeats(&repeats, seq_len, MIN_DISTANCE);
    println!("После фильтрации: {}", filtered.len());

    // --- Этап 2: карта рестрикции ---
    println!("\n--- Карта рестрикции ---");
    let enzymes = restr::parse_enzymes("data/enzymes.csv")?;
    println!("Загружено ферментов: {}", enzymes.len());

    let mut restriction_map = restr::build_map(seq, &enzymes);
    restriction_map.retain(|_, sites| sites.len() <= 500);

    // Нормализуем имена
    let mut normalized_map: HashMap<String, Vec<usize>> = HashMap::new();
    for (name, sites) in &restriction_map {
        let short = normalize_name(name);
        normalized_map.entry(short).or_default().extend(sites);
    }

    let mut with_sites: Vec<_> = normalized_map
        .iter()
        .filter(|(_, sites)| !sites.is_empty())
        .collect();
    with_sites.sort_by_key(|(_, sites)| -(sites.len() as i64));

    println!(
        "Ферментов (нормализовано, ≤500 сайтов): {}",
        with_sites.len()
    );
    for (name, sites) in with_sites.iter().take(10) {
        println!("  {}: {} сайтов", name, sites.len());
    }

    // --- Этап 3: диагностические фрагменты ---
    println!("\n--- Диагностические фрагменты ---");
    let mut all_solutions: Vec<(usize, String, String, char, usize)> = Vec::new();

    for (idx, repeat) in filtered.iter().enumerate() {
        let solutions = find_diagnostic_fragments(repeat, &normalized_map, seq_len);
        for sol in solutions {
            all_solutions.push((
                idx,
                sol.enzyme1,
                sol.enzyme2,
                sol.ring,
                sol.fragment_length,
            ));
        }
    }

    println!("Всего диагностических вариантов: {}", all_solutions.len());

    let mut unique: HashSet<(usize, String, String, char)> = HashSet::new();
    for (idx, e1, e2, ring, _) in &all_solutions {
        unique.insert((*idx, e1.clone(), e2.clone(), *ring));
    }
    println!(
        "Уникальных комбинаций (повтор + пара ферментов): {}",
        unique.len()
    );

    // Группировка по парам ферментов
    let mut enzyme_pair_coverage: HashMap<(String, String), HashSet<usize>> = HashMap::new();
    for (idx, e1, e2, _ring, _) in &all_solutions {
        let key = if e1 <= e2 {
            (e1.clone(), e2.clone())
        } else {
            (e2.clone(), e1.clone())
        };
        enzyme_pair_coverage.entry(key).or_default().insert(*idx);
    }

    // Сортируем по покрытию
    let mut ranked: Vec<_> = enzyme_pair_coverage
        .iter()
        .map(|((e1, e2), repeats)| (e1.clone(), e2.clone(), repeats.len()))
        .collect();
    ranked.sort_by_key(|(_, _, count)| -(*count as i64));

    println!("\n=== Топ-20 пар ферментов по покрытию ===");
    println!(
        "{:<20} | {:<20} | {}",
        "enzyme1", "enzyme2", "покрыто повторов"
    );
    println!("{:-<20}-+-{:-<20}-+-{:-<8}", "", "", "");
    for (e1, e2, count) in ranked.iter().take(20) {
        println!("{:<20} | {:<20} | {}", e1, e2, count);
    }

    // --- Этап 4: жадный отбор ---
    println!("\n--- Жадный отбор оптимального набора ---");
    let selected = greedy_select(&enzyme_pair_coverage, filtered.len());

    println!("Отобрано реакций: {}", selected.len());
    println!(
        "{:<3} | {:<20} | {:<20} | {:<8} | {}",
        "№", "фермент 1", "фермент 2", "покрыто", "накопл."
    );
    println!("{:-<3}-+-{:-<20}-+-{:-<20}-+-{:-<8}-+-{:-<8}", "", "", "", "", "");
    let mut cumulative: HashSet<usize> = HashSet::new();
    for (i, (e1, e2)) in selected.iter().enumerate() {
        if let Some(repeats) = enzyme_pair_coverage.get(&(e1.clone(), e2.clone())) {
            cumulative.extend(repeats);
            println!(
                "{:<3} | {:<20} | {:<20} | {:<8} | {}/{}",
                i + 1,
                e1,
                e2,
                repeats.len(),
                cumulative.len(),
                filtered.len()
            );
        }
    }

    // --- Экспорт в CSV ---
    println!("\n--- Экспорт ---");
    export_results(
        "diagnostic_results.csv",
        &all_solutions,
        &enzyme_pair_coverage,
        &selected,
        filtered.len(),
    )?;
    println!("Результаты сохранены в diagnostic_results.csv");

    println!("\nОбщее время: {:.2} сек", start.elapsed().as_secs_f64());
    Ok(())
}

// ---- Структуры ----

#[derive(Debug, Clone)]
struct Repeat {
    start1: usize,
    end1: usize,
    start2: usize,
    end2: usize,
    length: usize,
}

#[derive(Debug, Clone)]
struct DiagnosticSolution {
    start1: usize,
    end1: usize,
    start2: usize,
    end2: usize,
    ring: char,
    enzyme1: String,
    enzyme2: String,
    fragment_length: usize,
}

// ---- Вспомогательные функции ----

fn find_fasta(dir: &str) -> Result<std::path::PathBuf> {
    let dir = Path::new(dir);
    if !dir.exists() {
        anyhow::bail!("Папка '{}' не найдена.", dir.display());
    }
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if let Some(ext) = path.extension() {
            if ext == "fasta" || ext == "fa" || ext == "fna" {
                return Ok(path);
            }
        }
    }
    anyhow::bail!("В папке '{}' не найдено .fasta или .fa файлов.", dir.display())
}

fn hash_kmer(kmer: &[u8]) -> u64 {
    let mut h: u64 = 0;
    for &b in kmer.iter().take(16) {
        h = h.wrapping_mul(5).wrapping_add(match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        });
    }
    h
}

fn filter_repeats(repeats: &[Repeat], genome_len: usize, min_dist: usize) -> Vec<Repeat> {
    let mut sorted = repeats.to_vec();
    sorted.sort_by_key(|r| -(r.length as i64));

    let mut kept: Vec<Repeat> = Vec::new();
    for r in &sorted {
        let distance = distance_between(r, genome_len);
        if distance < min_dist {
            continue;
        }
        let overlaps = kept.iter().any(|existing| {
            overlaps_any(existing.start1, existing.end1, r.start1, r.end1)
                || overlaps_any(existing.start2, existing.end2, r.start2, r.end2)
        });
        if !overlaps {
            kept.push(r.clone());
        }
    }
    kept.sort_by_key(|r| r.start1);
    kept
}

fn distance_between(r: &Repeat, genome_len: usize) -> usize {
    let d1 = if r.start2 >= r.end1 {
        r.start2 - r.end1
    } else {
        genome_len - r.end1 + r.start2
    };
    let d2 = if r.start1 >= r.end2 {
        r.start1 - r.end2
    } else {
        genome_len - r.end2 + r.start1
    };
    d1.min(d2)
}

fn overlaps_any(start1: usize, end1: usize, start2: usize, end2: usize) -> bool {
    start1 < end2 && start2 < end1
}

fn find_diagnostic_fragments(
    repeat: &Repeat,
    rmap: &HashMap<String, Vec<usize>>,
    genome_len: usize,
) -> Vec<DiagnosticSolution> {
    let mut solutions = Vec::new();
    let flank = FLANK;

    let region_a1_start = repeat.start1.saturating_sub(flank);
    let region_a1_end = (repeat.start1 + flank).min(genome_len);
    let region_a2_start = repeat.end1.saturating_sub(flank);
    let region_a2_end = (repeat.end1 + flank).min(genome_len);
    let region_b1_start = repeat.start2.saturating_sub(flank);
    let region_b1_end = (repeat.start2 + flank).min(genome_len);
    let region_b2_start = repeat.end2.saturating_sub(flank);
    let region_b2_end = (repeat.end2 + flank).min(genome_len);

    for (enzyme1, sites1) in rmap {
        let sites_near_a1: Vec<usize> = sites1
            .iter()
            .filter(|&&p| p >= region_a1_start && p <= region_a1_end)
            .copied()
            .collect();
        let sites_near_a2: Vec<usize> = sites1
            .iter()
            .filter(|&&p| p >= region_a2_start && p <= region_a2_end)
            .copied()
            .collect();

        for (enzyme2, sites2) in rmap {
            let sites_near_b1: Vec<usize> = sites2
                .iter()
                .filter(|&&p| p >= region_b1_start && p <= region_b1_end)
                .copied()
                .collect();
            let sites_near_b2: Vec<usize> = sites2
                .iter()
                .filter(|&&p| p >= region_b2_start && p <= region_b2_end)
                .copied()
                .collect();

            for &s1 in &sites_near_a1 {
                for &s2 in &sites_near_a2 {
                    for &s3 in &sites_near_b1 {
                        for &s4 in &sites_near_b2 {
                            let len_no_recomb_1 = s2.saturating_sub(s1);
                            let len_no_recomb_2 = s4.saturating_sub(s3);
                            let len_recomb_1 = s4.saturating_sub(s1);
                            let len_recomb_2 = s2.saturating_sub(s3);

                            if (len_no_recomb_1 != len_recomb_1
                                || len_no_recomb_2 != len_recomb_2)
                                && (len_recomb_1 >= 50 || len_recomb_2 >= 50)
                            {
                                solutions.push(DiagnosticSolution {
                                    start1: repeat.start1,
                                    end1: repeat.end1,
                                    start2: repeat.start2,
                                    end2: repeat.end2,
                                    ring: 'A',
                                    enzyme1: enzyme1.clone(),
                                    enzyme2: enzyme2.clone(),
                                    fragment_length: len_recomb_1.max(len_recomb_2),
                                });
                            }
                        }
                    }
                }
            }
        }
    }
    solutions
}

fn normalize_name(name: &str) -> String {
    let clean = name
        .trim()
        .replace("WarmStart ", "")
        .replace("WarmStart", "");
    let clean = clean.replace(" §", "").replace("§", "");
    let clean = clean
        .trim_end_matches("-HFv2")
        .trim_end_matches("-HF")
        .trim_end_matches("-v2");
    clean.trim().to_string()
}

fn greedy_select(
    coverage: &HashMap<(String, String), HashSet<usize>>,
    total_repeats: usize,
) -> Vec<(String, String)> {
    let mut selected: Vec<(String, String)> = Vec::new();
    let mut covered: HashSet<usize> = HashSet::new();

    let mut pairs: Vec<(&(String, String), &HashSet<usize>)> = coverage.iter().collect();
    pairs.sort_by_key(|(_, repeats)| -(repeats.len() as i64));

    for _ in 0..pairs.len() {
        if covered.len() >= total_repeats {
            break;
        }
        let mut best: Option<((String, String), usize)> = None;
        for ((e1, e2), repeats) in &pairs {
            let new_coverage = repeats.difference(&covered).count();
            if new_coverage > 0 {
                match best {
                    Some((_, best_count)) if new_coverage > best_count => {
                        best = Some(((e1.clone(), e2.clone()), new_coverage));
                    }
                    None => {
                        best = Some(((e1.clone(), e2.clone()), new_coverage));
                    }
                    _ => {}
                }
            }
        }
        if let Some(((e1, e2), _)) = best {
            selected.push((e1.clone(), e2.clone()));
            if let Some(repeats) = coverage.get(&(e1.clone(), e2.clone())) {
                covered.extend(repeats);
            }
        } else {
            break;
        }
    }
    selected
}

fn export_results(
    path: &str,
    all: &[(usize, String, String, char, usize)],
    coverage: &HashMap<(String, String), HashSet<usize>>,
    selected: &[(String, String)],
    total_repeats: usize,
) -> Result<()> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "repeat_idx,enzyme1,enzyme2,ring,fragment_length,in_optimal_set"
    )?;

    let optimal_set: HashSet<(String, String)> = selected.iter().cloned().collect();

    for (idx, e1, e2, ring, frag_len) in all {
        let key = if e1 <= e2 {
            (e1.clone(), e2.clone())
        } else {
            (e2.clone(), e1.clone())
        };
        let in_opt = optimal_set.contains(&key);
        writeln!(
            file,
            "{},{},{},{},{},{}",
            idx, e1, e2, ring, frag_len, in_opt
        )?;
    }

    writeln!(file, "\n--- Optimal set summary ---")?;
    writeln!(file, "enzyme1,enzyme2,repeats_covered")?;
    let mut cumulative: HashSet<usize> = HashSet::new();
    for (e1, e2) in selected {
        if let Some(repeats) = coverage.get(&(e1.clone(), e2.clone())) {
            cumulative.extend(repeats);
            writeln!(file, "{},{},{}", e1, e2, repeats.len())?;
        }
    }
    writeln!(
        file,
        "\nTotal repeats covered: {}/{}",
        cumulative.len(),
        total_repeats
    )?;

    Ok(())
}
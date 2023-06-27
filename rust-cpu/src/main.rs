use chemfiles::Frame;
use chemfiles::Trajectory;
use indicatif::ProgressBar;
use std::fs::File;
use std::io::{self, BufWriter, Write};

fn calculate_distance_matrix(distance_matrix: &mut Vec<Vec<f64>>, frame: &Frame) {
    let n_atoms = distance_matrix.len();

    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            let dist = frame.distance(i, j);
            distance_matrix[i][j] = dist;
            distance_matrix[j][i] = dist;
        }
    }
}

fn distance_cutoff(distance: f64, rcutoff: f64) -> f64 {
    if distance > rcutoff {
        0.0
    } else {
        0.5 * (f64::cos(std::f64::consts::PI * distance / rcutoff) + 1.0)
    }
}

fn compute_g2_element(distance: f64, eta: f64, rcutoff: f64, rshift: f64) -> f64 {
    if distance > 0.0 {
        f64::exp(-eta * (distance - rshift).powi(2)) * distance_cutoff(distance, rcutoff)
    } else {
        0.0
    }
}

fn compute_g2(distances: &[f64], eta: f64, rcutoff: f64, rshift: f64) -> f64 {
    distances
        .iter()
        .map(|&d| compute_g2_element(d, eta, rcutoff, rshift))
        .sum()
}

fn build_g2_matrix(distance_matrix: &[Vec<f64>], eta: f64, rcutoff: f64, rshift: f64) -> Vec<f64> {
    distance_matrix
        .iter()
        .map(|distances| compute_g2(distances, eta, rcutoff, rshift))
        .collect()
}

fn compute_mean_g2(
    trajectory: &mut chemfiles::Trajectory,
    eta: f64,
    rcutoff: f64,
    rshift: f64,
    n_atoms: usize
) -> f64 {
    let mut frame = Frame::new();
    let mut distance_matrix = vec![vec![0.0; n_atoms]; n_atoms];

    let mut mean_g2_value = 0.0;

    let file_length = trajectory.nsteps();
    for _i in 0..file_length {
        trajectory.read(&mut frame).unwrap();
        calculate_distance_matrix(&mut distance_matrix, &frame);
        let g2_matrix = build_g2_matrix(&distance_matrix, eta, rcutoff, rshift);
        mean_g2_value += g2_matrix.iter().sum::<f64>() / g2_matrix.len() as f64;
    }

    mean_g2_value
}

fn main() -> io::Result<()> {
    let filename = "100CH3OH-CG-200.xtc";
    let mut trajectory = Trajectory::open(&filename, 'r').unwrap();
    let mut frame = Frame::new();
    trajectory.read(&mut frame);
    let n_atoms = frame.size();
    let rcutoff = 10.0;
    let eta_range: Vec<f64> = (0..=1).map(|x| x as f64).collect();
    let rshift_range: Vec<f64> = (0..=1).map(|x| x as f64).collect();
    let n_eta = eta_range.len();
    let n_rshift = rshift_range.len();
    let mut g2_values_matrix = vec![vec![0.0; n_rshift]; n_eta];

    println!("Iterations: {}", n_eta * n_rshift);

    let bar = ProgressBar::new((n_eta * n_rshift) as u64);

    for (i, eta) in eta_range.iter().enumerate() {
        for (j, rshift) in rshift_range.iter().enumerate() {
            let mut trajectory = Trajectory::open(&filename, 'r').unwrap();
            g2_values_matrix[i][j] = compute_mean_g2(&mut trajectory, *eta, rcutoff, *rshift, n_atoms);
            bar.inc(1);
        }
    }

    bar.finish();

    let output_file = "rust_g2_values.txt";
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);

    for row in g2_values_matrix {
        let line: Vec<String> = row.iter().map(|&x| format!("{:.e}", x)).collect();
        let line = line.join("    ");
        writeln!(writer, "{}", line)?;
    }

    Ok(())
}

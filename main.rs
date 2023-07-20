// extern crate rand;
extern crate rand;
extern crate stats;

use std::cmp::Ordering;
use std::fs::File;
use std::io::Write;
use std::iter;
use std::path::Path;

use rand::prelude::ThreadRng;
use rand::Rng;

// use std::rand;
// use std::rand::Rng;


struct Problem {
    x: Vec<f64>,
    y: f64,
    lows: Vec<f64>,
    highs: Vec<f64>,
    fitness: Box<dyn Fn(&Vec<f64>) -> f64>,
    criteria: Box<dyn Fn(f64, f64) -> bool>,
    factible: Box<dyn Fn(&Vec<f64>, &Vec<f64>, &Vec<f64>) -> bool>,

}

struct SolverOpt {
    max_iter: u64,
    sampling: u64,
}


fn perturb(solution: &Vec<f64>, rng: &mut ThreadRng) -> Vec<f64> {
    let mut new_solution = solution.clone();
    for x in new_solution.iter_mut() {
        let perturbation: f64 = rng.gen_range(-0.1, 0.1);
        println!("pert {:?} x {:?}",perturbation, x);
        *x += perturbation;
        println!("xafter {:?}", x);
        // if fitnessg
    }
    return new_solution;
}

fn initial_solution(lows: &Vec<f64>, highs: &Vec<f64>, n: u64) -> Vec<f64> {
    let mut rng = rand::thread_rng();

    return (0..n)
        .map(|i| rng.gen_range(10.0,
                               10.1))
        .collect::<Vec<_>>();
}

fn quadratic(lows: &Vec<f64>, highs: &Vec<f64>) -> Problem {
    let n = 1;

    let initial = initial_solution(lows, highs, n);
    let fitness = |x: &Vec<f64>| { x[0] * x[0] };

    let initial_fitness = fitness(&initial);
    return Problem {
        x: initial,
        y: initial_fitness,
        lows: lows.clone(),
        highs: highs.clone(),
        fitness: Box::new(fitness),
        criteria: Box::new(criteria),
        factible: Box::new(factible),
    };
}

fn factible(x: &Vec<f64>, lows: &Vec<f64>, highs: &Vec<f64>) -> bool {
    for i in 0..x.len() {
        if x[i] < lows[i] &&
            x[i] > highs[i] {
            return false;
        }
    }
    return true;
}

fn criteria(old_y: f64, new_y: f64) -> bool {
    return new_y < old_y;
}

fn ackleys(initial: Option<Vec<f64>>, lows: &Vec<f64>, highs: &Vec<f64>) -> Problem {
    let n = 2;

    let real_inital = match initial {
        None => initial_solution(lows, highs, n),
        Some(n) => n
    };

    let fitness = |x: &Vec<f64>| {
        return -(x[1] + 47.0) * f64::sin(f64::sqrt(f64::abs(x[0] / 2.0 +
            x[1] + 47.0))) - x[0] * f64::sin(f64::sqrt(f64::abs
            (x[0] - (x[1] + 47.0))));
    };


    let initial_fitness = fitness(&real_inital);
    return Problem {
        x: real_inital,
        y: initial_fitness,
        lows: lows.clone(),
        highs: highs.clone(),
        fitness: Box::new(fitness),
        criteria: Box::new(criteria),
        factible: Box::new(factible),
    };
}

fn descent(opts: SolverOpt, problem: &mut Problem) -> Vec<f64> {
    let mut rng = rand::thread_rng();

    let mut iter = 0;
    let mut progress = Vec::with_capacity(opts.max_iter as usize);

    while iter < opts.max_iter {
        // Perturb
        let mut new_xs: Vec<Vec<f64>> = vec![];
        let mut new_ys: Vec<f64> = vec![];

        for _ in 0..opts.sampling {
            let new_solution = perturb(&problem.x, &mut rng);
            new_ys.push(((*problem).fitness)(&new_solution));
            new_xs.push(new_solution);
        }

        let mut best_y: f64 = new_ys[0];
        let mut best_x: &Vec<f64> = &new_xs[0];
        println!("{:?} {:?}", new_ys[0], new_xs[0]);
        for i in new_xs[1..].iter().enumerate() {
            if ((*problem).factible)(i.1, &problem.lows, &problem.highs) {
                if ((*problem).criteria)(best_y, new_ys[i.0 +1]) {
                    best_y = new_ys[i.0 +1];
                    best_x = i.1;
                }
            }
        }
        iter += 1;
    }

    return progress;
}

fn run_n(initial: Option<Vec<f64>>, which_problem: &str, n: u32,
         lows: &Vec<f64>, highs: &Vec<f64>, sampling: u64) {
    let mut ys: Vec<f64> = vec![];
    let mut xs: Vec<Vec<f64>> = vec![];

    for _ in 0..n {
        let mut problem;
        if which_problem == "ack" {
            problem = ackleys(initial.clone(), lows, highs);
        } else {
            problem = quadratic(lows, highs)
        }
        descent(SolverOpt {
            max_iter: 50,
            sampling: sampling,
        }, &mut problem);
        ys.push(problem.y);
        xs.push(problem.x);
    }
    let index_besty = ys.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(Ordering::Equal)).unwrap();


    println!(" {:?} {:5} {:.5} & {:.5} & {:.5} & {:.5}",
            xs[index_besty.0],
            index_besty.1,
             ys.clone().iter().cloned().fold(0. / 0., f64::min),
             ys.clone().iter().cloned().fold(0. / 0., f64::max),
             stats::mean(ys.clone().into_iter()),
             stats::stddev(ys.clone().into_iter()));
}
//
// fn chart() {
//     let mut problem = ackleys(None, );
//     let progress = descent(SolverOpt {
//         max_iter: 50,
//         sampling: 1,
//     }, &mut problem);
//     {
//         let mut file = File::create(Path::new("ackleys_100.csv")).unwrap();
//
//         for i in progress.iter().enumerate() {
//             write!(file, "{},{}\n", i.0, i.1).unwrap();
//         }
//     }
// }

fn main() {

    // println!("SDHC Quadratic");
    // run_n(None, "quad", 30, &vec![-10.0, -10.0],
    //       &vec![10.0, 10.0], 100);
    println!("HC Quadratic");
    run_n(None, "quad", 1, &vec![-10.0, -10.0],
          &vec![10.0, 10.0], 1);
}

use line_drawing::Bresenham;
use minifb::{Key, Window, WindowOptions, Scale};

const X_RESOLUTION : usize = 600;
const Y_RESOLUTION : usize = 150;
const SCALING : Scale = Scale::X2;

const FOUR_NINTHS : f32 = 4.0 / 9.0;
const ONE_NINTH : f32 = 1.0 / 9.0;
const ONE_THIRTYSIXTH : f32 = 1.0 / 36.0;

const U0 : f32 = 0.04; // flow speed
const VISCOSITY : f32 = 0.004;
const STEPS_PER_FRAME : i32 = 20;
const FPS : i32 = 60;

const CONSTRAST : f32 = 5.0;
const N_COLORS : usize = 400;


// 0 1 2
// 3 4 5
// 6 7 8

fn main() {
    let mut time : i32 = 0;
    let mut step_count : i32 = 0;

    let mut n_nw = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_n = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_ne = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_w = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_o = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_e = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_sw = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_s = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut n_se = vec![0.0; X_RESOLUTION * Y_RESOLUTION];

    let mut rho = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut ux = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut uy = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut curl = vec![0.0; X_RESOLUTION * Y_RESOLUTION];
    let mut barrier = vec![false; X_RESOLUTION * Y_RESOLUTION];

    let barrier_size = 24;
    let range = (Y_RESOLUTION / 2 - barrier_size)..(Y_RESOLUTION / 2 + barrier_size + 1);
    for y in range {
        let x = (Y_RESOLUTION as f32 / 3.0).round() as usize;
        // barrier[x + y * X_RESOLUTION] = true;

        // barrier[x + 32 + (y - 8) * X_RESOLUTION] = true;
    }

    let xl = X_RESOLUTION as i32;
    let yh = Y_RESOLUTION as i32;
    let x1 = (Y_RESOLUTION as f32 / 3.0).round() as i32;
    for (x, y) in Bresenham::new((0, 0), (x1, yh * 2 / 3)) {
        barrier[x as usize + y as usize * X_RESOLUTION] = true;
        barrier[(x + 1) as usize + y as usize * X_RESOLUTION] = true;
    }

    for (x, y) in Bresenham::new((xl * 2 / 3, yh - 1), (xl * 3 / 4, yh * 1 / 4)) {
        barrier[x as usize + y as usize * X_RESOLUTION] = true;
        barrier[(x + 1) as usize + y as usize * X_RESOLUTION] = true;
    }

    for (x, y) in Bresenham::new((xl * 1 / 3, yh / 2), (xl * 1 / 2, yh / 2)) {
        barrier[x as usize + y as usize * X_RESOLUTION] = true;
        barrier[x as usize + (y + 1) as usize * X_RESOLUTION] = true;
    }



    let mut red_color_list = [0; N_COLORS + 2];
    let mut green_color_list = [0; N_COLORS + 2];
    let mut blue_color_list = [0; N_COLORS + 2];
    for i in 0..N_COLORS + 1 {
        let r;
        let g;
        let b;

        let float_index = i as f32;
        let eighth_n_colors = N_COLORS as f32 / 8.0;
        let fourth_n_colors = N_COLORS as f32 / 4.0;

        if float_index < eighth_n_colors {
            r = 0;
            g = 0;
            b = (255.0 * (float_index + eighth_n_colors) / fourth_n_colors).round() as i32;
        } else if float_index < 3.0 * eighth_n_colors {
            r = 0;
            g = (255.0 * (float_index - eighth_n_colors) / fourth_n_colors).round() as i32;
            b = 255;
        } else if float_index < 5.0 * eighth_n_colors {
            r = (255.0 * (float_index - 3.0 * eighth_n_colors) / fourth_n_colors).round() as i32;
            g = 255;
            b = 255 - r;
        } else if float_index < 7.0 * eighth_n_colors {
            r = 255;
            g = (255.0 * (-float_index + 7.0 * eighth_n_colors) / fourth_n_colors).round() as i32;
            b = 0;
        } else {
            r = (255.0 * (-float_index + 9.0 * eighth_n_colors) / fourth_n_colors).round() as i32;
            g = 0;
            b = 0;
        }

        red_color_list[i] = r;
        green_color_list[i] = g;
        blue_color_list[i] = b;

    }

    println!("Hello, world!");
    
    init_fluid(&mut curl, &mut rho, &mut ux, &mut uy, &mut n_nw, &mut n_n, &mut n_ne, &mut n_w, &mut n_o, &mut n_e, &mut n_sw, &mut n_s, &mut n_se);

    let mut window = init_canvas();

    while window.is_open() && !window.is_key_down(Key::Escape) {
        // println!("{:?}", nE);
        simulate(&mut rho, &mut ux, &mut uy, &mut barrier, &mut time, &mut step_count, &mut n_nw, &mut n_n, &mut n_ne, &mut n_w, &mut n_o, &mut n_e, &mut n_sw, &mut n_s, &mut n_se);
        paint_canvas(&mut curl, &mut window, &mut barrier, &mut ux, &mut uy, red_color_list, green_color_list, blue_color_list);
    }
    
}

fn init_canvas() -> Window {
    let mut window = Window::new(
        "Test - ESC to exit",
        X_RESOLUTION,
        Y_RESOLUTION,
        WindowOptions {
            resize: true,
            scale: SCALING,
            ..WindowOptions::default()
        }
    )
    .unwrap_or_else(|e| {
        panic!("{}", e);
    });

    window.limit_update_rate(Some(std::time::Duration::from_micros(1000000 / FPS as u64)));

    return window;
}

fn paint_canvas(
    curl : &mut Vec<f32>, window : &mut Window, barrier : &mut Vec<bool>, ux : &mut Vec<f32>, uy : &mut Vec<f32>,
    red_color_list : [i32; N_COLORS + 2], green_color_list : [i32; N_COLORS + 2], blue_color_list : [i32; N_COLORS + 2],

) {
    let mut buffer: Vec<u32> = vec![0; X_RESOLUTION * Y_RESOLUTION];
    let contrast = f32::powf(1.2, CONSTRAST);
    let mut c_index;

    compute_curl(curl, ux, uy);
                
    for y in 0..Y_RESOLUTION {
        for x in 0..X_RESOLUTION {
            if barrier[x + y * X_RESOLUTION] {
                c_index = N_COLORS + 1;
            } else {
                c_index = ((N_COLORS as f32) * (curl[x + y * X_RESOLUTION] * 5.0 * contrast + 0.5)).round() as usize;
                
                if c_index > N_COLORS {
                    c_index = N_COLORS;
                }
            }

            let red_component = (red_color_list[c_index] * 256 * 256) as u32;
            let green_component = (green_color_list[c_index] * 256) as u32;
            let blue_component = blue_color_list[c_index] as u32;

            buffer[x + y * X_RESOLUTION] = red_component + green_component + blue_component;
            // buffer[x + y * X_RESOLUTION] = ((red_color_list[c_index] << 16) | (green_color_list[c_index] << 8) | blue_color_list[c_index]) as u32;
        }
    }
    
    /*for i in 0..buffer.len(){
        buffer[i] = (i as u32 % 255) * 255 * 255 + 120 * 255 + 255;
    }*/

    // We unwrap here as we want this code to exit if it fails. Real applications may want to handle this in a different way
    window
        .update_with_buffer(&buffer, X_RESOLUTION, Y_RESOLUTION)
        .unwrap();
}

fn compute_curl(curl : &mut Vec<f32>, ux : &mut Vec<f32>, uy : &mut Vec<f32>) {
    for y in 1..Y_RESOLUTION-1 {
        for x in 1..X_RESOLUTION-1 {
            curl[x + y * X_RESOLUTION] = uy[x + 1 + y * X_RESOLUTION] - uy[x - 1 + y * X_RESOLUTION] - ux[x + (y + 1) * X_RESOLUTION] + ux[x + (y - 1) * X_RESOLUTION];
        }
    }
}

fn collide(
    rho : &mut Vec<f32>, ux : &mut Vec<f32>, uy : &mut Vec<f32>,
    n_nw : &mut Vec<f32>, n_n : &mut Vec<f32>, n_ne : &mut Vec<f32>, 
    n_w : &mut Vec<f32>, n_o : &mut Vec<f32>, n_e : &mut Vec<f32>, 
    n_sw : &mut Vec<f32>, n_s : &mut Vec<f32>, n_se : &mut Vec<f32>
) {
    let omega = 1.0 / (3.0 * VISCOSITY + 0.5);

    for y in 1..Y_RESOLUTION - 1 {
        for x in 1 ..X_RESOLUTION - 1 {
            let i = x + y * X_RESOLUTION;

            let this_rho = n_nw[i] + n_n[i] + n_ne[i] + n_w[i] + n_o[i] + n_e[i] + n_sw[i] + n_s[i] + n_se[i];
            rho[i] = this_rho;
            let this_ux = (n_e[i] + n_ne[i] + n_se[i] - n_w[i] - n_nw[i] - n_sw[i]) / this_rho;
            ux[i] = this_ux;
            let this_uy = (n_n[i] + n_nw[i] + n_ne[i] - n_s[i] - n_sw[i] - n_se[i]) / this_rho;
            uy[i] = this_uy;

            let one9thrho = ONE_NINTH * this_rho;		// pre-compute a bunch of stuff for optimization
            let one36thrho = ONE_THIRTYSIXTH * this_rho;
            let ux3 = 3.0 * this_ux;
            let uy3 = 3.0 * this_uy;
            let ux2 = this_ux * this_ux;
            let uy2 = this_uy * this_uy;
            let uxuy2 = 2.0 * this_ux * this_uy;
            let u2 = ux2 + uy2;
            let u215 = 1.5 * u2;

            n_o[i] += omega * (FOUR_NINTHS * this_rho * (1.0                    - u215) - n_o[i]);
            n_e[i] += omega * (   one9thrho * (1.0 + ux3       + 4.5*ux2        - u215) - n_e[i]);
            n_w[i] += omega * (   one9thrho * (1.0 - ux3       + 4.5*ux2        - u215) - n_w[i]);
            n_n[i] += omega * (   one9thrho * (1.0 + uy3       + 4.5*uy2        - u215) - n_n[i]);
            n_s[i] += omega * (   one9thrho * (1.0 - uy3       + 4.5*uy2        - u215) - n_s[i]);
            n_ne[i] += omega * (  one36thrho * (1.0 + ux3 + uy3 + 4.5*(u2+uxuy2) - u215) - n_ne[i]);
            n_se[i] += omega * (  one36thrho * (1.0 + ux3 - uy3 + 4.5*(u2-uxuy2) - u215) - n_se[i]);
            n_nw[i] += omega * (  one36thrho * (1.0 - ux3 + uy3 + 4.5*(u2-uxuy2) - u215) - n_nw[i]);
            n_sw[i] += omega * (  one36thrho * (1.0 - ux3 - uy3 + 4.5*(u2+uxuy2) - u215) - n_sw[i]);
        }
    }

    for y in 1..Y_RESOLUTION-2 {
        let to_coord = X_RESOLUTION - 1 + y * X_RESOLUTION;
        let from_coord = X_RESOLUTION - 2 + y * X_RESOLUTION;

        n_w[to_coord] = n_w[from_coord];
        n_nw[to_coord] = n_nw[from_coord];
        n_sw[to_coord] = n_sw[from_coord];
    }
}

fn stream(
    barrier : &mut Vec<bool>,
    n_nw : &mut Vec<f32>, n_n : &mut Vec<f32>, n_ne : &mut Vec<f32>, 
    n_w : &mut Vec<f32>, _n_o : &mut Vec<f32>, n_e : &mut Vec<f32>, 
    n_sw : &mut Vec<f32>, n_s : &mut Vec<f32>, n_se : &mut Vec<f32>
) {

    let north_calc = || {
        for y in (1..Y_RESOLUTION-1).rev() {
            for x in 1..X_RESOLUTION-1 {
                n_n[x + y * X_RESOLUTION] = n_n[x + (y - 1) * X_RESOLUTION];
                n_nw[x + y * X_RESOLUTION] = n_nw[x + 1 + (y - 1) * X_RESOLUTION];
            }
        }
    };

    let east_calc = || {
        for y in (1..Y_RESOLUTION-1).rev() {
            for x in (1..X_RESOLUTION-1).rev() {
                n_e[x + y * X_RESOLUTION] = n_e[x - 1 + y * X_RESOLUTION];
                n_ne[x + y * X_RESOLUTION] = n_ne[x - 1 + (y - 1) * X_RESOLUTION];
            }
        }
    };

    rayon::join(north_calc, east_calc);

    let south_calc = || {
        for y in 1..Y_RESOLUTION-1 {
            for x in (1..X_RESOLUTION-1).rev() {
                n_s[x + y * X_RESOLUTION] = n_s[x + (y + 1) * X_RESOLUTION];
                n_se[x + y * X_RESOLUTION] = n_se[x - 1 + (y + 1) * X_RESOLUTION];
            }
        }
    };

    let west_calc = || {
        for y in 1..Y_RESOLUTION-1 {
            for x in 1..X_RESOLUTION-1 {
                n_w[x + y * X_RESOLUTION] = n_w[x + 1 + y * X_RESOLUTION];
                n_sw[x + y * X_RESOLUTION] = n_sw[x + 1 + (y + 1) * X_RESOLUTION];
            }
        }
    };

    rayon::join(south_calc, west_calc);

    for y in 1..Y_RESOLUTION-1 {
        for x in 1..X_RESOLUTION-1 {
            if barrier[x + y * X_RESOLUTION] {
                let index = x +  y * X_RESOLUTION;

                n_e[x + 1 + y * X_RESOLUTION] = n_w[index];
                n_w[x - 1 + y * X_RESOLUTION] = n_e[index];
                n_n[x + (y + 1) * X_RESOLUTION] = n_s[index];
                n_s[x + (y - 1) * X_RESOLUTION] = n_n[index];

                n_ne[x + 1 + (y + 1) * X_RESOLUTION] = n_sw[index];
                n_nw[x - 1 + (y + 1) * X_RESOLUTION] = n_se[index];
                n_se[x + 1 + (y - 1) * X_RESOLUTION] = n_nw[index];
                n_sw[x - 1 + (y - 1) * X_RESOLUTION] = n_ne[index];
            }
        }
    }

}

fn simulate(
    rho : &mut Vec<f32>, ux : &mut Vec<f32>, uy : &mut Vec<f32>, barrier : &mut Vec<bool>,
    time : &mut i32, step_count : &mut i32,
    n_nw : &mut Vec<f32>, n_n : &mut Vec<f32>, n_ne : &mut Vec<f32>, 
    n_w : &mut Vec<f32>, n_o : &mut Vec<f32>, n_e : &mut Vec<f32>, 
    n_sw : &mut Vec<f32>, n_s : &mut Vec<f32>, n_se : &mut Vec<f32>
) {
    set_boundaries(rho, ux, uy, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);

    for _step in 0..STEPS_PER_FRAME {
        collide(rho, ux, uy, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);
        stream(barrier, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);
        *time += 1;
    }

    *step_count += STEPS_PER_FRAME;

    let mut stable = true;
    for x in 0..X_RESOLUTION {
        let index = x + (Y_RESOLUTION / 2) * X_RESOLUTION;
        if rho[index] <= 0.0 {
            stable = false;
        }
    }
    if !stable {
        println!("Unstable!");
    }
}

fn set_boundaries(
    rho : &mut Vec<f32>, ux : &mut Vec<f32>, uy : &mut Vec<f32>,
    n_nw : &mut Vec<f32>, n_n : &mut Vec<f32>, n_ne : &mut Vec<f32>, 
    n_w : &mut Vec<f32>, n_o : &mut Vec<f32>, n_e : &mut Vec<f32>, 
    n_sw : &mut Vec<f32>, n_s : &mut Vec<f32>, n_se : &mut Vec<f32>
){
    for x in 0..X_RESOLUTION {
        set_equil(x, 0, U0, 0.0, 1.0, rho, ux, uy, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);
        set_equil(x, Y_RESOLUTION - 1, U0, 0.0, 1.0, rho, ux, uy, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);
    }

    for y in 1..Y_RESOLUTION-1 {
        set_equil(0, y, U0, 0.0, 1.0, rho, ux, uy, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);
        set_equil(X_RESOLUTION - 1, y, U0, 0.0, 1.0, rho, ux, uy, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);
    }
}

fn init_fluid(
    curl : &mut Vec<f32>,
    rho : &mut Vec<f32>, ux : &mut Vec<f32>, uy : &mut Vec<f32>,
    n_nw : &mut Vec<f32>, n_n : &mut Vec<f32>, n_ne : &mut Vec<f32>, 
    n_w : &mut Vec<f32>, n_o : &mut Vec<f32>, n_e : &mut Vec<f32>, 
    n_sw : &mut Vec<f32>, n_s : &mut Vec<f32>, n_se : &mut Vec<f32>
){
    for y in 0..Y_RESOLUTION {
        for x in 0..X_RESOLUTION {
            set_equil(x, y, U0, 0.0, 1.0, rho, ux, uy, n_nw, n_n, n_ne, n_w, n_o, n_e, n_sw, n_s, n_se);
            curl[x + y * X_RESOLUTION] = 0.0;
        }
    }
}

fn set_equil(
    x: usize, y: usize, newux: f32, newuy: f32, newrho: f32, 
    rho : &mut Vec<f32>, ux : &mut Vec<f32>, uy : &mut Vec<f32>,
    n_nw : &mut Vec<f32>, n_n : &mut Vec<f32>, n_ne : &mut Vec<f32>, 
    n_w : &mut Vec<f32>, n_o : &mut Vec<f32>, n_e : &mut Vec<f32>, 
    n_sw : &mut Vec<f32>, n_s : &mut Vec<f32>, n_se : &mut Vec<f32>
) {
    let i = x + y * X_RESOLUTION;

    let ux3 = 3.0 * newux;
    let uy3 = 3.0 * newuy;
    let ux2 = newux * newux;
    let uy2 = newuy * newuy;
    let uxuy2 = 2.0 * newux * newuy;
    let u2 = ux2 + uy2;
    let u215 = 1.5 * u2;

    n_o[i] = FOUR_NINTHS * newrho * (1.0                              - u215);
    n_e[i] =   ONE_NINTH * newrho * (1.0 + ux3       + 4.5*ux2        - u215);
    n_w[i] =   ONE_NINTH * newrho * (1.0 - ux3       + 4.5*ux2        - u215);
    n_n[i] =   ONE_NINTH * newrho * (1.0 + uy3       + 4.5*uy2        - u215);
    n_s[i] =   ONE_NINTH * newrho * (1.0 - uy3       + 4.5*uy2        - u215);
    n_ne[i] =  ONE_THIRTYSIXTH * newrho * (1.0 + ux3 + uy3 + 4.5*(u2+uxuy2) - u215);
    n_se[i] =  ONE_THIRTYSIXTH * newrho * (1.0 + ux3 - uy3 + 4.5*(u2-uxuy2) - u215);
    n_nw[i] =  ONE_THIRTYSIXTH * newrho * (1.0 - ux3 + uy3 + 4.5*(u2-uxuy2) - u215);
    n_sw[i] =  ONE_THIRTYSIXTH * newrho * (1.0 - ux3 - uy3 + 4.5*(u2+uxuy2) - u215);

    rho[i] = newrho;
    ux[i] = newux;
    uy[i] = newuy;
}

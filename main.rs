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

struct Lattice {
    n_nw: Vec<f32>,
    n_n: Vec<f32>,
    n_ne: Vec<f32>,
    n_w: Vec<f32>,
    n_o: Vec<f32>,
    n_e: Vec<f32>,
    n_sw: Vec<f32>,
    n_s: Vec<f32>,
    n_se: Vec<f32>,

    rho: Vec<f32>,
    ux: Vec<f32>,
    uy: Vec<f32>,
    curl: Vec<f32>,
    barrier: Vec<bool>,
}

impl Lattice {
    fn init_fluid(&mut self) -> () {
        for y in 0..Y_RESOLUTION {
            for x in 0..X_RESOLUTION {
                self.set_equil(x, y, U0, 0.0, 1.0);
                self.curl[x + y * X_RESOLUTION] = 0.0;
            }
        }
    }

    fn set_equil(&mut self, x: usize, y: usize, newux: f32, newuy: f32, newrho: f32) -> () {
        let i = x + y * X_RESOLUTION;

        let ux3 = 3.0 * newux;
        let uy3 = 3.0 * newuy;
        let ux2 = newux * newux;
        let uy2 = newuy * newuy;
        let uxuy2 = 2.0 * newux * newuy;
        let u2 = ux2 + uy2;
        let u215 = 1.5 * u2;

        self.n_o[i] = FOUR_NINTHS * newrho * (1.0                              - u215);
        self.n_e[i] =   ONE_NINTH * newrho * (1.0 + ux3       + 4.5*ux2        - u215);
        self.n_w[i] =   ONE_NINTH * newrho * (1.0 - ux3       + 4.5*ux2        - u215);
        self.n_n[i] =   ONE_NINTH * newrho * (1.0 + uy3       + 4.5*uy2        - u215);
        self.n_s[i] =   ONE_NINTH * newrho * (1.0 - uy3       + 4.5*uy2        - u215);
        self.n_ne[i] =  ONE_THIRTYSIXTH * newrho * (1.0 + ux3 + uy3 + 4.5*(u2+uxuy2) - u215);
        self.n_se[i] =  ONE_THIRTYSIXTH * newrho * (1.0 + ux3 - uy3 + 4.5*(u2-uxuy2) - u215);
        self.n_nw[i] =  ONE_THIRTYSIXTH * newrho * (1.0 - ux3 + uy3 + 4.5*(u2-uxuy2) - u215);
        self.n_sw[i] =  ONE_THIRTYSIXTH * newrho * (1.0 - ux3 - uy3 + 4.5*(u2+uxuy2) - u215);

        self.rho[i] = newrho;
        self.ux[i] = newux;
        self.uy[i] = newuy;
    }

    fn set_boundaries(&mut self) -> () {
        for x in 0..X_RESOLUTION {
            self.set_equil(x, 0, U0, 0.0, 1.0);
            self.set_equil(x, Y_RESOLUTION - 1, U0, 0.0, 1.0);
        }
    
        for y in 1..Y_RESOLUTION-1 {
            self.set_equil(0, y, U0, 0.0, 1.0);
            self.set_equil(X_RESOLUTION - 1, y, U0, 0.0, 1.0);
        }
    }

    fn add_line_boundary(&mut self, start: (i32, i32), end: (i32, i32), thick_x: bool, thick_y: bool) -> () {
        for (x, y) in Bresenham::new(start, end) {
            self.barrier[x as usize + y as usize * X_RESOLUTION] = true;
            if thick_x { 
                self.barrier[x as usize + 1 + y as usize * X_RESOLUTION] = true;
            }
            if thick_y {
                self.barrier[x as usize + (y + 1) as usize * X_RESOLUTION] = true;
            }
        }
    }

    fn simulate(&mut self) -> () {
        self.set_boundaries();

        for _step in 0..STEPS_PER_FRAME {
            self.collide();
            self.stream();
        }
    }

    fn collide(&mut self) -> () {
        let omega = 1.0 / (3.0 * VISCOSITY + 0.5);

        for y in 1..Y_RESOLUTION - 1 {
            for x in 1 ..X_RESOLUTION - 1 {
                let i = x + y * X_RESOLUTION;

                let this_rho = self.n_nw[i] + self.n_n[i] + self.n_ne[i] + self.n_w[i] + self.n_o[i] + self.n_e[i] + self.n_sw[i] + self.n_s[i] + self.n_se[i];
                let this_ux = (self.n_e[i] + self.n_ne[i] + self.n_se[i] - self.n_w[i] - self.n_nw[i] - self.n_sw[i]) / this_rho;
                let this_uy = (self.n_n[i] + self.n_nw[i] + self.n_ne[i] - self.n_s[i] - self.n_sw[i] - self.n_se[i]) / this_rho;
                
                self.rho[i] = this_rho;
                self.ux[i] = this_ux;
                self.uy[i] = this_uy;

                let one9thrho = ONE_NINTH * this_rho;		// pre-compute a bunch of stuff for optimization
                let one36thrho = ONE_THIRTYSIXTH * this_rho;
                let ux3 = 3.0 * this_ux;
                let uy3 = 3.0 * this_uy;
                let ux2 = this_ux * this_ux;
                let uy2 = this_uy * this_uy;
                let uxuy2 = 2.0 * this_ux * this_uy;
                let u2 = ux2 + uy2;
                let u215 = 1.5 * u2;

                self.n_o[i] += omega * (FOUR_NINTHS * this_rho * (1.0                    - u215) - self.n_o[i]);
                self.n_e[i] += omega * (   one9thrho * (1.0 + ux3       + 4.5*ux2        - u215) - self.n_e[i]);
                self.n_w[i] += omega * (   one9thrho * (1.0 - ux3       + 4.5*ux2        - u215) - self.n_w[i]);
                self.n_n[i] += omega * (   one9thrho * (1.0 + uy3       + 4.5*uy2        - u215) - self.n_n[i]);
                self.n_s[i] += omega * (   one9thrho * (1.0 - uy3       + 4.5*uy2        - u215) - self.n_s[i]);
                self.n_ne[i] += omega * (  one36thrho * (1.0 + ux3 + uy3 + 4.5*(u2+uxuy2) - u215) - self.n_ne[i]);
                self.n_se[i] += omega * (  one36thrho * (1.0 + ux3 - uy3 + 4.5*(u2-uxuy2) - u215) - self.n_se[i]);
                self.n_nw[i] += omega * (  one36thrho * (1.0 - ux3 + uy3 + 4.5*(u2-uxuy2) - u215) - self.n_nw[i]);
                self.n_sw[i] += omega * (  one36thrho * (1.0 - ux3 - uy3 + 4.5*(u2+uxuy2) - u215) - self.n_sw[i]);
            }
        }

        for y in 1..Y_RESOLUTION-2 {
            let to_coord = X_RESOLUTION - 1 + y * X_RESOLUTION;
            let from_coord = X_RESOLUTION - 2 + y * X_RESOLUTION;

            self.n_w[to_coord] = self.n_w[from_coord];
            self.n_nw[to_coord] = self.n_nw[from_coord];
            self.n_sw[to_coord] = self.n_sw[from_coord];
        }
    }

    fn stream(&mut self) -> () {
        let north_calc = || {
            for y in (1..Y_RESOLUTION-1).rev() {
                for x in 1..X_RESOLUTION-1 {
                    self.n_n[x + y * X_RESOLUTION] = self.n_n[x + (y - 1) * X_RESOLUTION];
                    self.n_nw[x + y * X_RESOLUTION] = self.n_nw[x + 1 + (y - 1) * X_RESOLUTION];
                }
            }
        };
    
        let east_calc = || {
            for y in (1..Y_RESOLUTION-1).rev() {
                for x in (1..X_RESOLUTION-1).rev() {
                    self.n_e[x + y * X_RESOLUTION] = self.n_e[x - 1 + y * X_RESOLUTION];
                    self.n_ne[x + y * X_RESOLUTION] = self.n_ne[x - 1 + (y - 1) * X_RESOLUTION];
                }
            }
        };
    
        rayon::join(north_calc, east_calc);
    
        let south_calc = || {
            for y in 1..Y_RESOLUTION-1 {
                for x in (1..X_RESOLUTION-1).rev() {
                    self.n_s[x + y * X_RESOLUTION] = self.n_s[x + (y + 1) * X_RESOLUTION];
                    self.n_se[x + y * X_RESOLUTION] = self.n_se[x - 1 + (y + 1) * X_RESOLUTION];
                }
            }
        };
    
        let west_calc = || {
            for y in 1..Y_RESOLUTION-1 {
                for x in 1..X_RESOLUTION-1 {
                    self.n_w[x + y * X_RESOLUTION] = self.n_w[x + 1 + y * X_RESOLUTION];
                    self.n_sw[x + y * X_RESOLUTION] = self.n_sw[x + 1 + (y + 1) * X_RESOLUTION];
                }
            }
        };
    
        rayon::join(south_calc, west_calc);
    
        for y in 1..Y_RESOLUTION-1 {
            for x in 1..X_RESOLUTION-1 {
                if self.barrier[x + y * X_RESOLUTION] {
                    let index = x +  y * X_RESOLUTION;
    
                    self.n_e[x + 1 + y * X_RESOLUTION] = self.n_w[index];
                    self.n_w[x - 1 + y * X_RESOLUTION] = self.n_e[index];
                    self.n_n[x + (y + 1) * X_RESOLUTION] = self.n_s[index];
                    self.n_s[x + (y - 1) * X_RESOLUTION] = self.n_n[index];
    
                    self.n_ne[x + 1 + (y + 1) * X_RESOLUTION] = self.n_sw[index];
                    self.n_nw[x - 1 + (y + 1) * X_RESOLUTION] = self.n_se[index];
                    self.n_se[x + 1 + (y - 1) * X_RESOLUTION] = self.n_nw[index];
                    self.n_sw[x - 1 + (y - 1) * X_RESOLUTION] = self.n_ne[index];
                }
            }
        }
    }

    fn compute_curl(&mut self) -> () {
        for y in 1..Y_RESOLUTION-1 {
            for x in 1..X_RESOLUTION-1 {
                self.curl[x + y * X_RESOLUTION] = self.uy[x + 1 + y * X_RESOLUTION] - self.uy[x - 1 + y * X_RESOLUTION] - self.ux[x + (y + 1) * X_RESOLUTION] + self.ux[x + (y - 1) * X_RESOLUTION];
            }
        }
    }

    fn init_canvas(&mut self) -> Window {
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

    fn paint_canvas(&mut self, color_bank : &ColorBank, window : &mut Window) -> () {
        let mut buffer: Vec<u32> = vec![0; X_RESOLUTION * Y_RESOLUTION];
        let contrast = f32::powf(1.2, CONSTRAST);
        let mut c_index;
    
        self.compute_curl();
                    
        for y in 0..Y_RESOLUTION {
            for x in 0..X_RESOLUTION {
                if self.barrier[x + y * X_RESOLUTION] {
                    c_index = N_COLORS + 1;
                } else {
                    c_index = ((N_COLORS as f32) * (self.curl[x + y * X_RESOLUTION] * 5.0 * contrast + 0.5)).round() as usize;
                    
                    if c_index > N_COLORS {
                        c_index = N_COLORS;
                    }
                }
    
                let red_component = (color_bank.red_color_list[c_index] * 256 * 256) as u32;
                let green_component = (color_bank.green_color_list[c_index] * 256) as u32;
                let blue_component = color_bank.blue_color_list[c_index] as u32;
    
                buffer[x + y * X_RESOLUTION] = red_component + green_component + blue_component;
                // buffer[x + y * X_RESOLUTION] = ((red_color_list[c_index] << 16) | (green_color_list[c_index] << 8) | blue_color_list[c_index]) as u32;
            }
        }
        
        window
            .update_with_buffer(&buffer, X_RESOLUTION, Y_RESOLUTION)
            .unwrap(); 
    }
}

#[derive(Debug)]
struct ColorBank {
    red_color_list: [i32; N_COLORS + 2],
    green_color_list: [i32; N_COLORS + 2],
    blue_color_list: [i32; N_COLORS + 2],
}

impl ColorBank {
    fn init_colors(&mut self) -> () {
        for i in 0..N_COLORS+1 {
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

            self.red_color_list[i] = r;
            self.green_color_list[i] = g;
            self.blue_color_list[i] = b;
        }
    }
}

fn main() {
    let mut lattice = Lattice {
        n_nw: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_n: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_ne: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_w: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_o: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_e: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_sw: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_s: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        n_se: vec![0.0; X_RESOLUTION * Y_RESOLUTION],

        rho: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        ux: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        uy: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        curl: vec![0.0; X_RESOLUTION * Y_RESOLUTION],
        barrier: vec![false; X_RESOLUTION * Y_RESOLUTION],
    };

    let mut color_bank = ColorBank {
        red_color_list: [0; N_COLORS + 2],
        green_color_list: [0; N_COLORS + 2],
        blue_color_list: [0; N_COLORS + 2],
    };
    color_bank.init_colors();

    lattice.init_fluid();

    let xl = X_RESOLUTION as i32;
    let yh = Y_RESOLUTION as i32;
    let x1 = (Y_RESOLUTION as f32 / 3.0).round() as i32;
    lattice.add_line_boundary((0, 0), (x1, yh * 2 / 3), true, false);
    lattice.add_line_boundary((xl * 2 / 3, yh - 1), (xl * 3 / 4, yh * 1 / 4), true, false);
    lattice.add_line_boundary((xl * 1 / 3, yh / 2), (xl * 1 / 2, yh / 2), false, true);

    let mut window = lattice.init_canvas();
    while window.is_open() && !window.is_key_down(Key::Escape) {
        lattice.simulate();
        lattice.paint_canvas(&color_bank, &mut window);
    }
}

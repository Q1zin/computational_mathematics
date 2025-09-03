static MAX_COUNT_LOOP: i32 = 1000;
static ROUND_COUNT: u32 = 3;

// ax + b = 0
struct LinearDegree {
    a: f32,
    b: f32,
    roots: Option<f32>,
}

impl LinearDegree {
    pub fn new(a: f32, b: f32) -> Self {
        Self { a, b, roots: None }
    }

    pub fn get_roots(&mut self) -> f32 {
        if self.a == 0. {
            panic!("a должно быть не равно 0!");
        }

        self.roots = Some(round_to(-self.b / self.a, ROUND_COUNT));

        self.roots.unwrap()
    }
}

// ax^2 + bx + c = 0
struct SecondDegree {
    a: f32,
    b: f32,
    c: f32,
    eps: f32,
    discr: Option<f32>,
    roots: Option<Vec<f32>>,
}

impl SecondDegree {
    pub fn new(a: f32, b: f32, c: f32, eps: f32) -> Self {
        Self {
            a,
            b,
            c,
            eps,
            discr: None,
            roots: None,
        }
    }

    fn discriminant(&self) -> f32 {
        if self.a == 0. {
            panic!("a должно быть не равно 0!");
        }
        self.b * self.b - 4.0 * self.a * self.c
    }

    pub fn get_discr(&mut self) -> f32 {
        match self.discr {
            Some(discr) => discr,
            None => {
                self.discr = Some(self.discriminant());
                return self.discr.unwrap();
            }
        }
    }

    pub fn get_roots(&mut self) -> Option<Vec<f32>> {
        if self.get_discr() < 0. {
            self.roots = None;
            return self.roots.clone();
        }

        if self.discr.unwrap() < -self.eps {
            self.roots = None
        } else if self.discr.unwrap().abs() < self.eps {
            let x = round_to(-self.b / (2. * self.a), ROUND_COUNT);
            self.roots = Some(vec![x])
        } else {
            let discr_sqrt = self.discr.unwrap().sqrt();
            let x1 = round_to((-self.b + discr_sqrt) / (2.*self.a), ROUND_COUNT);
            let x2 = round_to((-self.b - discr_sqrt) / (2.*self.a), ROUND_COUNT);
            self.roots = Some(vec![x1, x2])
        }

        return self.roots.clone();
    }

    pub fn get_derivative(&self) -> LinearDegree {
        LinearDegree::new(2. * self.a, self.b)
    }

    pub fn get_value(&self, x: f32) -> f32 {
        self.a * x * x + self.b * x + self.c
    }
}

struct ThirdDegree {
    a: f32,
    b: f32,
    c: f32,
    eps: f32,
    delta: f32,
    roots: Option<Vec<f32>>,
}

impl ThirdDegree {
    pub fn new(a: f32, b: f32, c: f32, eps: f32, delta: f32) -> Self {
        Self { a, b, c, eps, delta, roots: None }
    }

    pub fn get_derivative(&self) -> SecondDegree {
        SecondDegree::new(3., 2. * self.a, self.b, self.eps)
    }

    pub fn get_roots(&mut self) -> Option<Vec<f32>> {
        let mut first_derivative = self.get_derivative();
        let discr = first_derivative.get_discr();

        println!("discr {discr}");

        if discr < -self.eps {
            let value_zero = first_derivative.get_value(0.);
            let left: f32;
            let right: f32;

            if value_zero < -self.eps {
                // идём от 0 до +inf
                (left, right) = self.get_local_root(0., f32::INFINITY);
            } else if value_zero > self.eps {
                // идём от -inf до 0
                (left, right) = self.get_local_root(f32::NEG_INFINITY, 0.);
            } else {
                // то корень = 0
                (left, right) = (0., 0.)
            }

            println!("left {left} right {right}");

            let x = self.found_root(left, right);

            if self.get_derivative().get_value(x).abs() > self.eps {
                println!("Нашли корень {x} и ещё существует 2 комплексных корня");
            } else {
                println!("Нашли корень {x}, его кратность 3");
            }
        } else if discr >= self.eps {
            let critical_points = first_derivative.get_roots();
            match critical_points {             
                Some(v) => {
                    for x in v.iter() {
                        println!("{x} ");
                    }
                },
                None => {
                    print!("нету корней");
                }
            }
        } else {
            panic!("Прекол, хз что делать =)))")
        }

        return self.roots.clone();
    }

    fn get_local_root(&self, a: f32, b: f32) -> (f32, f32) {
        let go_left_or_right;

        if a == f32::NEG_INFINITY {
            go_left_or_right = -1;
        } else if b == f32::INFINITY {
            go_left_or_right = 1;
        } else {
            panic!("Ни какая из крайних точек не inf")
        }

        let mut left: f32;
        let mut right: f32;

        (left, right) = match go_left_or_right {
            -1 => (b - self.delta, b),
            1 => (a, a + self.delta),
            _ => panic!("Неверное значение go_left_or_right")
        };

        let mut count_loop = 0;
       
        while (self.get_value(left) * self.get_value(right)) > 0. && count_loop < MAX_COUNT_LOOP {
            count_loop += 1;
            println!("{count_loop} {left} {right}");
            (left, right) = match go_left_or_right {
                -1 => (left - self.delta, left),
                1 => (right, right + self.delta),
                _ => panic!("Неверное значение go_left_or_right")
            };
        }

        if count_loop == MAX_COUNT_LOOP {
            panic!("Не удалось дойти до корня(((");
        }

        return (left, right);
    }

    fn found_root(&self, left: f32, right: f32) -> f32 {
        let mut left = left;
        let mut right = right;

        let count_loop = 0;

        loop {
            let value_left = self.get_value(left);
            let value_right = self.get_value(right);

            // println!("value_left {value_left} value_right {value_right} left {left} right {right}");

            if value_left.abs() < self.eps && value_right.abs() < self.eps {
                return (left + right) / 2.;
            }

            if value_left.abs() <= self.eps {
                return left;
            }

            if value_right.abs() <= self.eps {
                return right;
            }

            let center: f32 = (left + right) / 2.;
            let value_center = self.get_value(center);

            if value_center.abs() < self.eps {
                return center;
            }

            if value_left * value_center < 0. {
                right = center;
            } else {
                left = center;
            }

            if count_loop >= MAX_COUNT_LOOP {
                panic!("Не удалось дойти до корня на отрезке")
            }
        }
    }

    pub fn get_value(&self, x: f32) -> f32 {
        x * x * x + self.a * x * x + self.b * x + self.c
    }
}

fn round_to(x: f32, decimals: u32) -> f32 {
    let factor = 10_f32.powi(decimals as i32);
    (x * factor).round() / factor
}

fn main() {
    let mut gg = ThirdDegree::new(-6., 12., -8., 0., 0.3);
    gg.get_roots();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roots_linear() {
        let mut equation_1 = LinearDegree::new(2.0, 10.0);
        assert_eq!(equation_1.get_roots(), -5.);

        let mut equation_2 = LinearDegree::new(-10.0, 2.0);
        assert_eq!(equation_2.get_roots(), 0.2);

        let mut equation_3 = LinearDegree::new(16.0, -4.0);
        assert_eq!(equation_3.get_roots(), 0.25);

        let mut equation_4 = LinearDegree::new(3.0, 10.0);
        assert_eq!(equation_4.get_roots(), -3.333);
    }

    #[test]
    fn test_discriminant_second_degree() {
        let eps = 1e-6;

        // D = b^2 - 4ac = 4 - 8 = -4
        let mut equation_1 = SecondDegree::new(2.0, 2.0, 1.0, eps);
        assert_eq!(equation_1.get_discr(), -4.);

        // D = (-3)^2 - 4*1*2 = 9 - 8 = 1
        let mut equation_2 = SecondDegree::new(1.0, -3.0, 2.0, eps);
        assert_eq!(equation_2.get_discr(), 1.);

        // D = 2^2 - 4*1*1 = 4 - 4 = 0
        let mut equation_3 = SecondDegree::new(1.0, 2.0, 1.0, eps);
        assert_eq!(equation_3.get_discr(), 0.);

        // D = 0^2 - 4*1*1 = -4
        let mut equation_4 = SecondDegree::new(1.0, 0.0, 1.0, eps);
        assert_eq!(equation_4.get_discr(), -4.);

        // D = (-7)^2 - 4*3.5*1 = 49 - 14 = 35
        let mut equation_5 = SecondDegree::new(3.5, -7.0, 1.0, eps);
        assert_eq!(equation_5.get_discr(), 35.);
    }

    #[test]
    fn test_roots_second_degree() {
        let eps = 1e-6;

        let mut equation_1 = SecondDegree::new(2.0, 2.0, 1.0, eps);
        assert_eq!(equation_1.get_roots(), None);

        let mut equation_2 = SecondDegree::new(1.0, -3.0, 2.0, eps);
        assert_eq!(equation_2.get_roots(), Some(vec![2.0, 1.0]));

        let mut equation_3 = SecondDegree::new(1.0, 2.0, 1.0, eps);
        assert_eq!(equation_3.get_roots(), Some(vec![-1.0]));

        let mut equation_4 = SecondDegree::new(1.0, 0.0, 1.0, eps);
        assert_eq!(equation_4.get_roots(), None);

        let mut equation_5 = SecondDegree::new(3.5, -7.0, 1.0, eps);
        assert_eq!(equation_5.get_roots(), Some(vec![1.845, 0.155]));
    }
}
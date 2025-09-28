use std::result;

use log::{info, error};
use env_logger::Env;

static MAX_COUNT_LOOP: i32 = 10000;
static ROUND_COUNT: u32 = 3;

// ax + b = 0
struct LinearEquation {
    a: f32,
    b: f32,
    roots: Option<f32>,
}

impl LinearEquation {
    pub fn new(a: f32, b: f32) -> Self {
        Self { a, b, roots: None }
    }

    pub fn get_value(&self, x: f32) -> f32 {
        x * self.a + self.b
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
struct QuadraticEquation {
    a: f32,
    b: f32,
    c: f32,
    eps: f32,
    discr: Option<f32>,
    roots: Option<Vec<f32>>,
}

impl QuadraticEquation {
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

    pub fn get_derivative(&self) -> LinearEquation {
        LinearEquation::new(2. * self.a, self.b)
    }

    pub fn get_value(&self, x: f32) -> f32 {
        self.a * x * x + self.b * x + self.c
    }
}

struct CubicEquation {
    a: f32,
    b: f32,
    c: f32,
    eps: f32,
    delta: f32,
    roots: Option<Vec<f32>>,
}

impl CubicEquation {
    pub fn new(a: f32, b: f32, c: f32, eps: f32, delta: f32) -> Self {
        Self { a, b, c, eps, delta, roots: None }
    }

    pub fn get_derivative(&self) -> QuadraticEquation {
        QuadraticEquation::new(3., 2. * self.a, self.b, self.eps)
    }

    pub fn get_roots(&mut self) -> Option<Vec<f32>> {
        info!("get_roots start");
        let mut first_derivative = self.get_derivative();
        let discr = first_derivative.get_discr();

        info!("discr {discr}");

        if discr <= self.eps {
            let value_zero = self.get_value(0.);
            info!("value_zero {value_zero}");
            let left: f32;
            let right: f32;

            if value_zero < -self.eps {
                (left, right) = self.get_local_root(0., f32::INFINITY);
            } else if value_zero > self.eps {
                (left, right) = self.get_local_root(f32::NEG_INFINITY, 0.);
            } else {
                (left, right) = (0., 0.)
            }

            info!("left {left} right {right}");

            let x = round_to(self.found_root(left, right), ROUND_COUNT);
            
            self.roots = Some(vec![x]);
        } else {
            let Some(mut critical_points) = first_derivative.get_roots() else {
                panic!("Нету корней, хотя должны быть");
            };

            critical_points.sort_by(|a, b| a.partial_cmp(b).unwrap());
            info!("critical_points {:?}", critical_points);

            if critical_points.len() != 2 {
                info!("Количество критических точек не равно 2, значит это монотонная функция =)))");
                critical_points.append(&mut vec![critical_points[0]]);
            }

            let mut left = critical_points[0];
            let mut right = critical_points[1];
            let value_left = self.get_value(left);
            let value_right = self.get_value(right);

            info!("value_left {value_left} right_left {value_right}");

            if value_left < -self.eps && value_right < -self.eps {
                info!("1");
                (left, right) = self.get_local_root(right, f32::INFINITY);
                let x = round_to(self.found_root(left, right), ROUND_COUNT);
                self.roots = Some(vec![x]);
            } else if value_left > self.eps && value_right > self.eps {
                info!("2");
                (left, right) = self.get_local_root(f32::NEG_INFINITY, left);
                let x = round_to(self.found_root(left, right), ROUND_COUNT);
                self.roots = Some(vec![x]);
            } else if value_left > self.eps && value_right.abs() < self.eps {
                info!("3");
                let (left, right_new) = self.get_local_root(f32::NEG_INFINITY, left);
                let x = round_to(self.found_root(left, right_new), ROUND_COUNT);
                self.roots = Some(vec![x, right]);
            } else if value_left.abs() < self.eps && value_right < -self.eps {
                info!("4");
                let (left_new, right) = self.get_local_root(right, f32::INFINITY);
                let x = round_to(self.found_root(left_new, right), ROUND_COUNT);
                self.roots = Some(vec![left, x]);
            } else if value_left > self.eps && value_right < -self.eps {
                info!("5 {left} {right}");
                let (left_1, right_1) = self.get_local_root(f32::NEG_INFINITY, left);
                let x1 = round_to(self.found_root(left_1, right_1), ROUND_COUNT);
                let x2 = round_to(self.found_root(left, right), ROUND_COUNT);
                let (left_3, right_3) = self.get_local_root(right, f32::INFINITY);
                let x3 = round_to(self.found_root(left_3, right_3), ROUND_COUNT);
                self.roots = Some(vec![x1, x2, x3]);
            } else if value_left.abs() < self.eps && value_right.abs() < self.eps{
                info!("6");
                self.roots = Some(vec![round_to((left + right) / 2., ROUND_COUNT)])
            } else {
                error!("Не все ситуации учел в случая расположения корней");
            }
        }

        info!("get_roots end");
        return self.roots.clone();
    }

    pub fn get_stat_roots(&mut self) -> Option<Vec<(f32, i32)>> {
        let Some(roots) = self.get_roots() else {
            return None;
        };

        let mut result: Vec<(f32, i32)> = vec![];

        for root in roots {
            let first_val = self.get_derivative().get_value(root);
            let second_val = self.get_derivative().get_derivative().get_value(root);

            let mut count = 1;
            if first_val.abs() <= self.eps {
                count += 1;
                if second_val.abs() <= self.eps {
                    count += 1;
                }
            }
            info!("root: {root}, count: {count}");

            result.push((root, count));
        }

        if result.len() == 0 {
            return None;
        }

        return Some(result);
    }

    fn get_local_root(&self, a: f32, b: f32) -> (f32, f32) {
        info!("get_local_root start");
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

        info!("0 {left} {right}");

        let mut count_loop = 0;
       
        while (self.get_value(left) * self.get_value(right)) > 0. && count_loop < MAX_COUNT_LOOP {
            count_loop += 1;
            info!("{count_loop} {left} {right}");
            (left, right) = match go_left_or_right {
                -1 => (left - self.delta, left),
                1 => (right, right + self.delta),
                _ => panic!("Неверное значение go_left_or_right")
            };
            info!("{count_loop} {left} {right}");
        }

        if count_loop == MAX_COUNT_LOOP {
            panic!("Не удалось дойти до корня(((");
        }

        info!("get_local_root end");

        return (left, right);
    }

    fn found_root(&self, left: f32, right: f32) -> f32 {
        info!("found_root start");
        let mut left = left;
        let mut right = right;

        info!("0 {left} {right}");

        let mut count_loop = 0;

        loop {
            count_loop += 1;
            let value_left = self.get_value(left);
            let value_right = self.get_value(right);

            info!("value_left {value_left} value_right {value_right} left {left} right {right}");

            if value_left.abs() < self.eps && value_right.abs() < self.eps {
                info!("found_root end 1 {}", (left + right) / 2.);
                return (left + right) / 2.;
            }

            if value_left.abs() <= self.eps {
                info!("found_root end 2 {}", left);
                return left;
            }

            if value_right.abs() <= self.eps {
                info!("found_root end 3 {}", right);
                return right;
            }

            let center: f32 = (left + right) / 2.;
            let value_center = self.get_value(center);
            if value_center.abs() < self.eps {
                info!("found_root end 4 {}", center);
                return center;
            }

            if value_left * value_center < 0. {
                right = center;
            } else {
                left = center;
            }

            info!("{count_loop} [{left} {right}] f(left)={value_left} f(right)={value_right} f(center)={value_center}");

            if count_loop >= MAX_COUNT_LOOP {
                panic!("Не удалось дойти до корня на отрезке с данной точностью")
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
    env_logger::Builder::from_env(Env::default().default_filter_or("error")).format_timestamp_millis().init();

    // let mut test = CubicEquation::new(-6., 12., -8., 1e-4, 0.3);
    let mut test = CubicEquation::new(3., 3., 1., 1e-5, 1.);
    println!("Итоговый ответ: {:?}", test.get_stat_roots());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roots_linear() {
        let mut equation_1 = LinearEquation::new(2.0, 10.0);
        assert_eq!(equation_1.get_roots(), -5.);

        let mut equation_2 = LinearEquation::new(-10.0, 2.0);
        assert_eq!(equation_2.get_roots(), 0.2);

        let mut equation_3 = LinearEquation::new(16.0, -4.0);
        assert_eq!(equation_3.get_roots(), 0.25);

        let mut equation_4 = LinearEquation::new(3.0, 10.0);
        assert_eq!(equation_4.get_roots(), -3.333);
    }

    #[test]
    fn test_discriminant_second_degree() {
        let eps = 1e-6;

        // D = b^2 - 4ac = 4 - 8 = -4
        let mut equation_1 = QuadraticEquation::new(2.0, 2.0, 1.0, eps);
        assert_eq!(equation_1.get_discr(), -4.);

        // D = (-3)^2 - 4*1*2 = 9 - 8 = 1
        let mut equation_2 = QuadraticEquation::new(1.0, -3.0, 2.0, eps);
        assert_eq!(equation_2.get_discr(), 1.);

        // D = 2^2 - 4*1*1 = 4 - 4 = 0
        let mut equation_3 = QuadraticEquation::new(1.0, 2.0, 1.0, eps);
        assert_eq!(equation_3.get_discr(), 0.);

        // D = 0^2 - 4*1*1 = -4
        let mut equation_4 = QuadraticEquation::new(1.0, 0.0, 1.0, eps);
        assert_eq!(equation_4.get_discr(), -4.);

        // D = (-7)^2 - 4*3.5*1 = 49 - 14 = 35
        let mut equation_5 = QuadraticEquation::new(3.5, -7.0, 1.0, eps);
        assert_eq!(equation_5.get_discr(), 35.);
    }

    #[test]
    fn test_roots_second_degree() {
        let eps = 1e-6;

        let mut equation_1 = QuadraticEquation::new(2.0, 2.0, 1.0, eps);
        assert_eq!(equation_1.get_roots(), None);

        let mut equation_2 = QuadraticEquation::new(1.0, -3.0, 2.0, eps);
        assert_eq!(equation_2.get_roots(), Some(vec![2.0, 1.0]));

        let mut equation_3 = QuadraticEquation::new(1.0, 2.0, 1.0, eps);
        assert_eq!(equation_3.get_roots(), Some(vec![-1.0]));

        let mut equation_4 = QuadraticEquation::new(1.0, 0.0, 1.0, eps);
        assert_eq!(equation_4.get_roots(), None);

        let mut equation_5 = QuadraticEquation::new(3.5, -7.0, 1.0, eps);
        assert_eq!(equation_5.get_roots(), Some(vec![1.845, 0.155]));
    }
}

#[cfg(test)]
mod third_degree_tests {
    use super::*;

    const EPS: f32 = 1e-4;
    const DELTA: f32 = 0.3;

    fn approx_eq(a: f32, b: f32, tol: f32) -> bool {
        (a - b).abs() <= tol
    }

    fn sort(v: &mut Vec<f32>) {
        v.sort_by(|x, y| x.partial_cmp(y).unwrap());
    }

    fn assert_vec_approx_eq_unordered(mut got: Vec<f32>, mut expected: Vec<f32>, tol: f32) {
        sort(&mut got);
        sort(&mut expected);
        assert_eq!(got.len(), expected.len(), "different lengths: got={got:?}, expected={expected:?}");
        for (g, e) in got.iter().zip(expected.iter()) {
            assert!(approx_eq(*g, *e, tol), "element mismatch: got={g}, expected={e}");
        }
    }

    // 1) discr < 0 и f(0) < 0  → корень единственный существует
    #[test]
    fn td_monotone_right_interval() {
        // x^3 + 2x - 1  (a=0, b=2, c=-1) → производная 3x^2 + 2a x + b = 3x^2 + 2
        let mut f = CubicEquation::new(0.0, 2.0, -1.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 1, "expected single root, got {:?}", got);
        assert!(approx_eq(got[0], 0.453, 1e-3), "root mismatch: got={}, expected≈0.453", got[0]);
    }

    // 2) discr < 0 и f(0) > 0  → корень единственный существует
    #[test]
    fn td_monotone_left_interval() {
        // x^3 + 2x + 1 → также монотонна
        let mut f = CubicEquation::new(0.0, 2.0, 1.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 1, "expected single root, got {:?}", got);
        assert!(approx_eq(got[0], -0.453, 1e-3), "root mismatch: got={}, expected≈-0.453", got[0]);
    }

    // 3) Ветка "1": value_left < 0 и value_right < 0 → один корень справа от правой крит. точки
    #[test]
    fn td_both_negative_one_root_right() {
        // x^3 - 4x - 6  → корень примерно 2.525102
        let mut f = CubicEquation::new(0.0, -4.0, -6.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 1, "expected single root, got {:?}", got);
        assert!(approx_eq(got[0], 2.5251, 1e-3), "root mismatch: got={}, expected≈2.5251", got[0]);
    }

    // 4) Ветка "2": value_left > 0 и value_right > 0 → один корень слева от левой крит. точки
    #[test]
    fn td_both_positive_one_root_left() {
        // x^3 - 4x + 6  → корень примерно -2.525102
        let mut f = CubicEquation::new(0.0, -4.0, 6.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 1, "expected single root, got {:?}", got);
        assert!(approx_eq(got[0], -2.5251, 1e-3), "root mismatch: got={}, expected≈-2.5251", got[0]);
    }

    // 5) Ветка "3": value_left > 0 и value_right ≈ 0 → правый крит. пункт — корень (двукратный),
    //    плюс ещё один корень слева
    #[test]
    fn td_right_critical_on_axis() {
        // x(x-2)^2 = x^3 - 4x^2 + 4x  → корни: 0 и 2 (2 — кратный)
        let mut f = CubicEquation::new(-4.0, 4.0, 0.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        // реализация возвращает [x_левый, right_critical]
        assert_vec_approx_eq_unordered(got, vec![0.0, 2.0], 2e-3);
    }

    // 6) Ветка "4": value_left ≈ 0 и value_right < 0 → левый крит. пункт — корень (двукратный),
    //    плюс ещё один корень справа
    #[test]
    fn td_left_critical_on_axis() {
        // (x-1)^2 (x-3) = x^3 - 5x^2 + 7x - 3 → корни: 1 (кратный), 3
        let mut f = CubicEquation::new(-5.0, 7.0, -3.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_vec_approx_eq_unordered(got, vec![1.0, 3.0], 2e-3);
    }

    // 7) Ветка "5": value_left > 0 и value_right < 0 → три действительных корня
    #[test]
    fn td_three_real_roots_symmetric() {
        // x(x^2 - 4) = x^3 - 4x → корни: -2, 0, 2
        let mut f = CubicEquation::new(0.0, -4.0, 0.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 3, "expected three roots, got {:?}", got);
        assert_vec_approx_eq_unordered(got, vec![-2.0, 0.0, 2.0], 2e-3);
    }

    // 8) Ещё пример ветки "5": три разных действительных корня без симметрии
    #[test]
    fn td_three_real_roots_asymmetric() {
        // (x+1)(x-2)(x-3) = x^3 - 4x^2 + x + 6 → корни: -1, 2, 3
        let mut f = CubicEquation::new(-4.0, 1.0, 6.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 3, "expected three roots, got {:?}", got);
        assert_vec_approx_eq_unordered(got, vec![-1.0, 2.0, 3.0], 2e-3);
    }

    // 9) Случай discr ≈ 0 (тройной корень): сейчас попадает в ветку error! и возвращает None
    #[test]
    fn td_triple_root_currently_none() {
        // (x-1)^3 = x^3 - 3x^2 + 3x - 1 → тройной корень 1
        let mut f = CubicEquation::new(-3.0, 3.0, -1.0, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 1, "expected three roots, got {:?}", got);
        assert_vec_approx_eq_unordered(got, vec![0.975], 2e-3); // но вообще говоря, если понизить DELTA, то будет 1.0
    }

    // 10) Три близких корня (проверяем устойчивость к точности/шагу delta)
    #[test]
    fn td_three_close_roots_precision() {
        // (x+2)(x-1)(x-1.1) = x^3 - 0.1x^2 - 3.1x + 2.2 → корни: -2, 1, 1.1
        let mut f = CubicEquation::new(-0.1, -3.1, 2.2, EPS, DELTA);
        let got = f.get_roots().expect("expected Some roots");
        assert_eq!(got.len(), 3, "expected three roots, got {:?}", got);
        assert_vec_approx_eq_unordered(got, vec![-2.0, 1.0, 1.1], 5e-3);
    }
}
use nalgebra::Vector2;

pub fn round_to_n_decimal_places(value: f64, decimals: u32) -> f64 {
    let factor = 10f64.powi(decimals as i32);
    (value * factor).round() / factor
}

pub fn approx_eq(a: &Vector2<f64>, b: &Vector2<f64>) -> bool {
    round_to_n_decimal_places(a[0], 6) == round_to_n_decimal_places(b[0], 6)
        && round_to_n_decimal_places(a[1], 6) == round_to_n_decimal_places(b[1], 6)
}

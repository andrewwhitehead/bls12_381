//! Document me

/// Evaluate isogeny map from curve with non-zero j-invariant.
pub trait IsogenyMap {
    /// Evaluate isogeny map.
    fn isogeny_map(&mut self);
}

// /// Trait executing Optimized Simplified SWU maps.
// ///
// /// Implemented for G1 and G2 based on https://eprint.iacr.org/2019/403.
// pub trait OsswuMap: CurveProjective {
//     /// Evaluate optimized simplified SWU map on supplied base field element.
//     fn osswu_map(u: &<Self as CurveProjective>::Base) -> Self;
// }

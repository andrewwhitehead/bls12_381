//! This module implements hash_to_curve, hash_to_field and related
//! hashing primitives for use with BLS signatures.

use core::ops::AddAssign;

pub(crate) mod chain;

mod expand_msg;
pub use self::expand_msg::{ExpandMessage, ExpandMsgXmd, ExpandMsgXof, InitExpandMessage};

mod g1;
mod g2;

use crate::generic_array::{typenum::Unsigned, ArrayLength, GenericArray};

/// Methods used in converting the output of hashed or encoded input
/// into an element of the field
pub trait MessageToField: for<'a> AddAssign<&'a Self> {
    /// The length of the data used to produce a component element
    type InputLength: ArrayLength<u8>;
    /// The component element
    type Pt: Copy + Default + core::fmt::Debug;

    /// Convert output keying material to a field element
    fn input_okm(okm: &GenericArray<u8, Self::InputLength>) -> Self::Pt;

    fn osswu_map(pt: &Self::Pt) -> Self;

    fn isogeny_map(&mut self);

    fn clear_h(&mut self);
}

/// Random oracle and injective maps to curve.
pub trait HashToCurve<'x, X>
where
    X: InitExpandMessage<'x>,
{
    /// Random oracle
    fn hash_to_curve(msg: impl AsRef<[u8]>, dst: &'x [u8]) -> Self;

    /// Injective encoding
    fn encode_to_curve(msg: impl AsRef<[u8]>, dst: &'x [u8]) -> Self;
}

impl<'x, F, X> HashToCurve<'x, X> for F
where
    F: MessageToField + core::fmt::Debug,
    X: InitExpandMessage<'x>,
{
    fn hash_to_curve(message: impl AsRef<[u8]>, dst: &'x [u8]) -> F {
        let mut p = {
            let mut u = [F::Pt::default(); 2];
            hash_to_field::<F, X>(message.as_ref(), dst, &mut u);
            // note: draft 10 performs isogeny_map for each component,
            // draft 7 performs one after adding the two. but adding two components
            // which aren't on the curve does not produce the expected result.
            let mut tmp = F::osswu_map(&u[0]);
            tmp.isogeny_map();
            let mut t2 = F::osswu_map(&u[1]);
            t2.isogeny_map();
            tmp.add_assign(&t2);
            tmp
        };
        p.clear_h();
        p
    }

    fn encode_to_curve(message: impl AsRef<[u8]>, dst: &'x [u8]) -> F {
        let mut u = [F::Pt::default(); 1];
        hash_to_field::<F, X>(message.as_ref(), dst, &mut u);
        println!("enc: {:?}", u);
        let mut p = F::osswu_map(&u[0]);
        println!("encb: {:?}", p);
        p.isogeny_map();
        println!("encc: {:?}", p);
        p.clear_h();
        println!("encd: {:?}", p);
        println!("encd: {:?}", p);
        p
    }
}

/// Hash to field for the type `F` using ExpandMessage variant `X`.
pub fn hash_to_field<'x, F, X>(
    message: &[u8],
    dst: &'x [u8],
    output: &mut [<F as MessageToField>::Pt],
) where
    F: MessageToField,
    X: InitExpandMessage<'x>,
{
    let len_per_elm = <F as MessageToField>::InputLength::to_usize();
    let len_in_bytes = output.len() * len_per_elm;
    let mut expander = X::init_expand(message, dst, len_in_bytes);

    let mut buf = GenericArray::<u8, <F as MessageToField>::InputLength>::default();
    for idx in 0..output.len() {
        expander.read_into(&mut buf[..]);
        println!("read buf: {}", hex::encode(&buf[..]));
        output[idx] = F::input_okm(&buf);
    }
}

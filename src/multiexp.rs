use core::ops;

use crypto_bigint::Limb;
use group::Group;
use subtle::{Choice, ConditionallySelectable};

use crate::scalar::{Inner as InnerScalar, Scalar};

/// Calculate a linear combination of points multiplied by scalars.
#[inline]
pub fn sum_of_products<P, Q>(mut points: &[P], mut scalars: &[Scalar]) -> Q
where
    P: Copy,
    Q: ConditionallySelectable + Group + for<'a> ops::Add<&'a P, Output = Q>,
{
    const BATCH: usize = 32;
    let mut remain = points.len();
    assert_eq!(remain, scalars.len());

    #[cfg(feature = "alloc")]
    if remain >= 32 {
        let mut norm = alloc::vec::Vec::with_capacity(remain);
        for i in 0..remain {
            norm.push(scalars[i].to_canonical());
        }
        return sum_of_products_batch(&points, &norm);
    }

    let mut tot = Q::identity();
    let mut norm = [InnerScalar::ZERO; BATCH];

    while remain > 0 {
        let c = remain.min(BATCH);
        for i in 0..c {
            norm[i] = scalars[i].to_canonical();
        }
        tot += sum_of_products_batch::<P, Q>(&points[..c], &norm[..c]);
        remain = remain - c;
        points = &points[c..];
        scalars = &scalars[c..];
    }

    tot
}

#[inline]
fn sum_of_products_batch<P, Q>(points: &[P], norm: &[InnerScalar]) -> Q
where
    P: Copy,
    Q: ConditionallySelectable + Group + for<'a> ops::Add<&'a P, Output = Q>,
{
    let mut acc = Q::identity();
    let mut i = Scalar::BIT_SIZE - 1;

    loop {
        for (j, point) in points.iter().enumerate() {
            acc.conditional_assign(&(acc + point), Choice::from(norm[j].bit_vartime(i) as u8));
        }
        if i == 0 {
            break;
        }
        i -= 1;
        acc = acc.double();
    }
    acc
}

/// Calculate a linear combination of points multiplied by scalars.
#[cfg(feature = "threaded")]
#[inline]
pub fn sum_of_products_vartime<P, Q>(points: &[P], scalars: &[Scalar]) -> Q
where
    P: Copy + Sync,
    Q: Group + From<P> + for<'a> ops::AddAssign<&'a P>,
{
    const BATCH: usize = 48;
    let remain = points.len();
    assert_eq!(remain, scalars.len());

    if remain > BATCH {
        let (l, r): (Q, Q) = rayon::join(
            || sum_of_products_vartime(&points[..BATCH], &scalars[..BATCH]),
            || sum_of_products_vartime(&points[BATCH..], &scalars[BATCH..]),
        );
        l + r
    } else {
        let mut norm = [Scalar::zero(); BATCH];
        for i in 0..remain {
            norm[i] = Scalar(scalars[i].to_canonical());
        }
        sum_of_products_pippenger(&points[..remain], &norm[..remain])
    }
}

/// Calculate a linear combination of points multiplied by scalars.
#[cfg(not(feature = "threaded"))]
#[inline]
pub fn sum_of_products_vartime<P, Q>(mut points: &[P], mut scalars: &[Scalar]) -> Q
where
    P: Copy + Sync,
    Q: Group + From<P> + for<'a> ops::AddAssign<&'a P> + for<'a> ops::AddAssign<&'a Q>,
{
    const BATCH: usize = 48;
    let mut remain = points.len();
    assert_eq!(remain, scalars.len());

    #[cfg(feature = "alloc")]
    if remain > BATCH {
        let mut norm = alloc::vec::Vec::with_capacity(remain);
        for i in 0..remain {
            norm.push(Scalar(scalars[i].to_canonical()));
        }
        return sum_of_products_pippenger(&points, &norm);
    }

    let mut tot = Q::identity();
    let mut norm = [Scalar::zero(); BATCH];

    while remain > 0 {
        let c = remain.min(BATCH);
        for i in 0..c {
            norm[i] = Scalar(scalars[i].to_canonical());
        }
        tot += sum_of_products_pippenger::<P, Q>(&points[..c], &norm[..c]);
        remain = remain - c;
        points = &points[c..];
        scalars = &scalars[c..];
    }

    tot
}

pub fn sum_of_products_pippenger<P, Q>(points: &[P], scalars: &[Scalar]) -> Q
where
    P: Copy,
    Q: Group + From<P> + for<'a> ops::AddAssign<&'a P>,
{
    use crypto_bigint::Word;

    const WINDOW: usize = 4;
    const NUM_BUCKETS: usize = (1 << WINDOW) - 1;
    const EDGE: usize = WINDOW - 1;

    // WINDOW is required to divide Word::BITS to avoid combining words
    debug_assert!((Word::BITS as usize) % WINDOW == 0);

    let num_components = points.len().min(scalars.len());
    let mut buckets = [Q::identity(); NUM_BUCKETS];
    let mut bit_sequence_index = InnerScalar::LIMBS * Limb::BIT_SIZE - 1;

    loop {
        let mut max_bucket = 1;
        let word_index = bit_sequence_index / (Word::BITS as usize);
        let bit_index = bit_sequence_index % (Word::BITS as usize);

        let shift = bit_index - EDGE;
        for i in 0..num_components {
            let bucket_index =
                ((scalars[i].0.as_words()[word_index] as usize) >> shift) & NUM_BUCKETS;
            if bucket_index > max_bucket {
                max_bucket = bucket_index;
                buckets[bucket_index - 1] = points[i].into();
            } else if bucket_index > 0 {
                buckets[bucket_index - 1] += &points[i];
            }
        }

        for i in (1..max_bucket).rev() {
            if i % 2 == 1 {
                buckets[i / 2] += buckets[i].double();
            } else {
                buckets[0] += buckets[i];
                buckets[i - 1] += buckets[i];
            }
            buckets[i] = Q::identity();
        }
        if bit_sequence_index < WINDOW {
            break buckets[0];
        }
        bit_sequence_index -= WINDOW;
        for _ in 0..WINDOW {
            buckets[0] = buckets[0].double();
        }
    }
}

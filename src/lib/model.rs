//! The two models a run scores calls under, and the distances they score by.
//!
//! The `chaff` model, the default, learns each filter's artifact fraction per
//! sample (see [`crate::prior`]) and scores end repair fill-in and copied
//! damage with a decay, `w(d) = exp(-d / s)`, from the end their artifact
//! favors, whose scale `s` it learns per library too: a polymerase fills
//! overhangs and copies lesions over lengths that vary from fragment to
//! fragment, so the evidence fades with distance rather than stopping at one.
//! End repair fill-in is strand-aware: a fill-in error sits on the strand the
//! polymerase extended, near that strand's 3' end. A-tailing keeps its
//! window, since its artifact is a non-templated A on the last base or two of
//! a 3' end.
//!
//! The `fgbio` model reproduces fgbio's `FilterSomaticVcf`: its per-call
//! `(2 * maf)^2` mutation prior and its windows from the nearest template end
//! for end repair fill-in and A-tailing. Copied damage, which fgbio lacks,
//! keeps its decay under fgbio's prior.

use std::fmt;
use std::str::FromStr;

use clap::ValueEnum;

/// The model that scores calls: its prior and the shape of its distance
/// models.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum Model {
    /// A learned prior, and learned decays for end repair fill-in and copied
    /// damage.
    #[default]
    Chaff,
    /// fgbio's per-call prior and windows, for parity with fgbio.
    Fgbio,
}

impl fmt::Display for Model {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Model::Chaff => write!(f, "chaff"),
            Model::Fgbio => write!(f, "fgbio"),
        }
    }
}

/// Where copied damage takes each call's artifact prior from under the
/// `chaff` model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum CopiedDamagePrior {
    /// On a BAM with single-strand consensus, the share of the library's
    /// positions as deep, with as many duplex changes, that chance explains
    /// (see [`crate::simplex`]); the learned fraction otherwise.
    #[default]
    Chance,
    /// The fraction learned per sample and stratum from the calls, even on a
    /// BAM with single-strand consensus, whose profile still fills the
    /// metrics.
    Learned,
}

impl fmt::Display for CopiedDamagePrior {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CopiedDamagePrior::Chance => write!(f, "chance"),
            CopiedDamagePrior::Learned => write!(f, "learned"),
        }
    }
}

/// A filter's distance: learned from the calls, or fixed in bases.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Distance {
    /// Learned per library, with the artifact fraction.
    #[default]
    Learned,
    /// A fixed number of bases.
    Bases(f64),
}

impl Distance {
    /// The fixed bases, or `learned` when there are none.
    pub fn bases(self, learned: f64) -> f64 {
        match self {
            Distance::Learned => learned,
            Distance::Bases(bases) => bases,
        }
    }
}

impl fmt::Display for Distance {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Distance::Learned => write!(f, "learned"),
            Distance::Bases(bases) => write!(f, "{bases}"),
        }
    }
}

impl FromStr for Distance {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        if text == "learned" {
            return Ok(Distance::Learned);
        }
        match text.parse::<f64>() {
            Ok(bases) if bases.is_finite() && bases > 0.0 => Ok(Distance::Bases(bases)),
            _ => Err(format!(
                "expected 'learned' or a positive number of bases, found: {text}"
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_a_distance_is_learned_or_a_positive_number_of_bases() {
        assert_eq!("learned".parse::<Distance>(), Ok(Distance::Learned));
        assert_eq!("12.5".parse::<Distance>(), Ok(Distance::Bases(12.5)));
        for text in ["0", "-3", "inf", "NaN", "Learned", ""] {
            assert!(text.parse::<Distance>().is_err(), "{text}");
        }
        assert_eq!(Distance::Learned.bases(21.0), 21.0);
        assert_eq!(Distance::Bases(15.0).bases(21.0), 15.0);
        assert_eq!(Distance::Bases(15.0).to_string(), "15");
    }
}

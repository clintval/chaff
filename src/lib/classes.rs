//! Substitution classes, lesion strands, and sequence context.
//!
//! A damage class such as `C>T` names the base on the damaged strand and what
//! that base reads as after polymerase copies it. A call explains a class on
//! either strand: `C>T` on the forward strand (REF `C`, ALT `T`) puts the lesion
//! on the forward strand, and `G>A` (its reverse complement) puts it on the
//! reverse strand.

use std::fmt;
use std::str::FromStr;

use anyhow::{bail, Result};

/// The complement of an upper-case DNA base; other bytes map to `N`.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Whether `base` is one of the four upper-case DNA bases.
fn is_dna(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

/// One strand of a double-stranded molecule, named by the reference strand it
/// matches.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Strand {
    /// The strand whose sequence matches the reference; its 5' end is the
    /// molecule's leftmost base.
    Forward,
    /// The complementary strand; its 5' end is the molecule's rightmost base.
    Reverse,
}

/// A damage class: a base on the damaged strand and the base it reads as.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct DamageClass {
    /// The undamaged base on the lesion strand, e.g. `C` for deamination.
    pub lesion: u8,
    /// The base the lesion reads as, e.g. `T` for deaminated cytosine.
    pub reads_as: u8,
}

impl DamageClass {
    /// Cytosine deamination, or deamination of 5-methylcytosine: `C>T`.
    pub const DEAMINATION: Self = Self {
        lesion: b'C',
        reads_as: b'T',
    };

    /// Guanine oxidation to 8-oxoguanine, which pairs with adenine: `G>T`.
    pub const OXIDATION: Self = Self {
        lesion: b'G',
        reads_as: b'T',
    };

    /// The strand carrying the lesion when this class explains the forward
    /// strand change `ref_base>alt_base`, or `None` when it does not.
    pub fn lesion_strand(&self, ref_base: u8, alt_base: u8) -> Option<Strand> {
        if ref_base == self.lesion && alt_base == self.reads_as {
            Some(Strand::Forward)
        } else if ref_base == complement(self.lesion) && alt_base == complement(self.reads_as) {
            Some(Strand::Reverse)
        } else {
            None
        }
    }

    /// The same class written from the other strand, e.g. `G>A` for `C>T`.
    pub fn reverse_complement(&self) -> Self {
        Self {
            lesion: complement(self.lesion),
            reads_as: complement(self.reads_as),
        }
    }
}

impl fmt::Display for DamageClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}>{}", self.lesion as char, self.reads_as as char)
    }
}

impl FromStr for DamageClass {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self> {
        let upper = s.trim().to_ascii_uppercase();
        let bytes = upper.as_bytes();
        if bytes.len() != 3 || bytes[1] != b'>' || !is_dna(bytes[0]) || !is_dna(bytes[2]) {
            bail!("a damage class is two bases joined by '>', quoted in a shell as 'C>T', but found: {s}");
        }
        if bytes[0] == bytes[2] {
            bail!("a damage class must change the base, but found: {s}");
        }
        Ok(Self {
            lesion: bytes[0],
            reads_as: bytes[2],
        })
    }
}

/// Reject a class list with a duplicate, counting a class and its reverse
/// complement as the same class.
pub fn validate_classes(classes: &[DamageClass]) -> Result<()> {
    for (i, a) in classes.iter().enumerate() {
        for b in &classes[i + 1..] {
            if a == b {
                bail!("damage class {a} is given twice");
            }
            if *a == b.reverse_complement() {
                bail!("damage classes {a} and {b} describe the same change on opposite strands");
            }
        }
    }
    Ok(())
}

/// Whether the lesion base sits in a CpG dinucleotide.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Context {
    /// The base is the C or the G of a CpG on the reference.
    CpG,
    /// Any other context, including an unknown neighbor.
    NonCpG,
}

impl Context {
    /// Classify the reference base at a position from its forward-strand
    /// neighbors. CpG is palindromic, so the result holds on both strands.
    pub fn of(prev: Option<u8>, base: u8, next: Option<u8>) -> Self {
        let base = base.to_ascii_uppercase();
        let prev = prev.map(|b| b.to_ascii_uppercase());
        let next = next.map(|b| b.to_ascii_uppercase());
        if (base == b'C' && next == Some(b'G')) || (base == b'G' && prev == Some(b'C')) {
            Context::CpG
        } else {
            Context::NonCpG
        }
    }
}

impl fmt::Display for Context {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Context::CpG => write!(f, "CpG"),
            Context::NonCpG => write!(f, "non-CpG"),
        }
    }
}

/// The pyrimidine-centric single-base substitution label (one of the six SBS
/// classes) for a forward-strand change, e.g. `C>T` for both `C>T` and `G>A`.
pub fn sbs6(ref_base: u8, alt_base: u8) -> Option<String> {
    let (r, a) = (ref_base.to_ascii_uppercase(), alt_base.to_ascii_uppercase());
    if !is_dna(r) || !is_dna(a) || r == a {
        return None;
    }
    let (r, a) = if matches!(r, b'C' | b'T') {
        (r, a)
    } else {
        (complement(r), complement(a))
    };
    Some(format!("{}>{}", r as char, a as char))
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case("C>T", b'C', b'T')]
    #[case("g>t", b'G', b'T')]
    #[case(" A>G ", b'A', b'G')]
    fn test_damage_class_parses(#[case] text: &str, #[case] lesion: u8, #[case] reads_as: u8) {
        let class: DamageClass = text.parse().unwrap();
        assert_eq!(class, DamageClass { lesion, reads_as });
    }

    #[rstest]
    #[case("C>C")]
    #[case("CT")]
    #[case("C>N")]
    #[case("CC>T")]
    fn test_damage_class_rejects_malformed(#[case] text: &str) {
        assert!(text.parse::<DamageClass>().is_err());
    }

    /// An unquoted `C>T` reaches the tool as `C`, the shell having taken `>T`
    /// as a redirection, so the error says to quote it.
    #[test]
    fn test_a_damage_class_cut_short_by_the_shell_says_to_quote_it() {
        let error = "C".parse::<DamageClass>().unwrap_err().to_string();
        assert!(error.contains("quoted in a shell as 'C>T'"), "{error}");
    }

    #[test]
    fn test_damage_class_display_round_trips() {
        assert_eq!(DamageClass::DEAMINATION.to_string(), "C>T");
        assert_eq!(DamageClass::OXIDATION.to_string(), "G>T");
    }

    #[rstest]
    #[case(DamageClass::DEAMINATION, b'C', b'T', Some(Strand::Forward))]
    #[case(DamageClass::DEAMINATION, b'G', b'A', Some(Strand::Reverse))]
    #[case(DamageClass::DEAMINATION, b'C', b'A', None)]
    #[case(DamageClass::OXIDATION, b'G', b'T', Some(Strand::Forward))]
    #[case(DamageClass::OXIDATION, b'C', b'A', Some(Strand::Reverse))]
    #[case(DamageClass::OXIDATION, b'G', b'A', None)]
    fn test_lesion_strand(
        #[case] class: DamageClass,
        #[case] ref_base: u8,
        #[case] alt_base: u8,
        #[case] expected: Option<Strand>,
    ) {
        assert_eq!(class.lesion_strand(ref_base, alt_base), expected);
    }

    #[test]
    fn test_validate_classes_rejects_reverse_complement_duplicates() {
        let classes = vec![DamageClass::DEAMINATION, "G>A".parse().unwrap()];
        let error = validate_classes(&classes).unwrap_err().to_string();
        assert!(
            error.contains("C>T and G>A describe the same change"),
            "{error}"
        );
        let twice = [DamageClass::DEAMINATION, DamageClass::DEAMINATION];
        let error = validate_classes(&twice).unwrap_err().to_string();
        assert_eq!(error, "damage class C>T is given twice");
        assert!(validate_classes(&[DamageClass::DEAMINATION, DamageClass::OXIDATION]).is_ok());
    }

    #[rstest]
    #[case(Some(b'A'), b'C', Some(b'G'), Context::CpG)]
    #[case(Some(b'C'), b'G', Some(b'A'), Context::CpG)]
    #[case(Some(b'a'), b'c', Some(b'g'), Context::CpG)]
    #[case(Some(b'G'), b'C', Some(b'A'), Context::NonCpG)]
    #[case(Some(b'T'), b'G', Some(b'C'), Context::NonCpG)]
    #[case(None, b'G', None, Context::NonCpG)]
    fn test_context(
        #[case] prev: Option<u8>,
        #[case] base: u8,
        #[case] next: Option<u8>,
        #[case] expected: Context,
    ) {
        assert_eq!(Context::of(prev, base, next), expected);
    }

    #[rstest]
    #[case(b'C', b'T', "C>T")]
    #[case(b'G', b'A', "C>T")]
    #[case(b'G', b'T', "C>A")]
    #[case(b'A', b'T', "T>A")]
    #[case(b'T', b'G', "T>G")]
    #[case(b'A', b'C', "T>G")]
    fn test_sbs6(#[case] ref_base: u8, #[case] alt_base: u8, #[case] expected: &str) {
        assert_eq!(sbs6(ref_base, alt_base).as_deref(), Some(expected));
    }

    #[test]
    fn test_sbs6_rejects_non_substitutions() {
        assert_eq!(sbs6(b'C', b'C'), None);
        assert_eq!(sbs6(b'N', b'C'), None);
    }
}

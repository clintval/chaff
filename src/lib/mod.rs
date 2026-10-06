//! `chaff`: separate somatic variant calls from library-preparation damage
//! artifacts in Duplex Sequencing and other UMI sequencing.
//!
//! # Differences from fgbio
//!
//! chaff ports the filters, likelihoods, and tests of fgbio's
//! `FilterSomaticVcf`, and with `--model fgbio` it writes fgbio 4.1.1's values
//! and FILTERs wherever overlapping mates agree in base and quality, no read
//! has an indel or soft clip between the call and its mate's 5' end, and every
//! base is Q2 or better. It matches fgbio where fgbio's choices are arbitrary:
//! a deletion at the site counts in the depth of its prior, a spanning
//! deletion `*` is no called allele, and an A-tailing site equally far from
//! both template ends is nearer the kept read's own end.
//!
//! It differs on purpose here:
//!
//! - **The chaff model, the default.** fgbio's mutation prior is near zero at
//!   duplex allele fractions, so alternate molecules inside the window make a
//!   call an artifact however many reference molecules sit there too; the
//!   chaff model learns the prior per library instead (see [`prior`]). fgbio's
//!   end repair window counts either template end and stops at a fixed
//!   distance; the chaff model measures from the 3' end of the strand each
//!   template was copied from and decays with distance at a learned scale (see
//!   [`model`]).
//! - **End repair.** A polymerase extends a recessed 3' end across a 5'
//!   overhang; fgbio's docs describe filling a 3' overhang.
//! - **Template ends.** Distances from both template ends count template
//!   bases: chaff walks both reads' CIGARs, the mate's from its `MC` tag, as
//!   [fgbio #1172](https://github.com/fulcrumgenomics/fgbio/pull/1172) does
//!   for clipping, so an indel counts by its length. Soft clips count and hard
//!   clips do not. fgbio measures the far end by insert size, and chaff never
//!   reads `TLEN`.
//! - **Overlapping mates.** They are called into one base: mates that agree
//!   keep the higher quality, and mates that disagree count as neither allele.
//!   fgbio keeps the first read of each name, so values differ where
//!   overlapping mates differ in base or quality.
//! - **Base errors.** A base's error probability is capped at 0.75, a random
//!   base's, so a Q0 or Q1 base cannot zero a likelihood.
//! - **Missing evidence.** A call with alternate but no reference molecules
//!   gets no INFO value; fgbio writes `NaN`.
//! - **Streaming.** The BAM is always streamed, never queried by index.
//! - **Number format.** Values keep htsjdk's rounding but are written in
//!   decimal: `0.00003218` for fgbio's `3.218e-05`.
#![warn(missing_docs)]

pub mod call;
pub mod classes;
pub mod copied_damage;
pub mod evidence;
pub mod filter;
pub mod io;
pub mod metrics;
pub mod model;
pub mod prior;
pub mod read_end;
pub mod reference;
pub mod testing;

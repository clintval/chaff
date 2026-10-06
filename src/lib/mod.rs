//! `chaff`: separate somatic variant calls from library-preparation damage
//! artifacts in duplex and UMI sequencing.
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

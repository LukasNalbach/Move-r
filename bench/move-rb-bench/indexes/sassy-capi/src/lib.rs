//! C interface to sassy's `Searcher`, used by move-rb-bench-apm.
//!
//! A searcher is created once per (profile, trace, all) combination and reused for all patterns
//! and all chunks of the text. The text is searched in overlapping chunks (overlap m + k - 1);
//! `sassy_capi_search` searches one chunk and returns only the number of matches.

use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use sassy::profiles::{Ascii, Dna, Iupac};
use sassy::Searcher;

/// A searcher for one of the three profiles.
pub enum SearcherType {
    Ascii(Searcher<Ascii>),
    Dna(Searcher<Dna>),
    Iupac(Searcher<Iupac>),
}

/// A searcher and its reporting mode.
pub struct Handle {
    searcher: SearcherType,
    /// true <=> all matches are reported, false <=> only the rightmost match of each local minimum
    all: bool,
}

/// Creates a searcher.
///
/// `profile` is "ascii", "dna" or "iupac". "dna" maps every byte to {A,C,G,T} (`(c >> 1) & 3`),
/// so a text containing a separator symbol has to be searched with "ascii".
/// `trace` selects whether the matches are computed with a traceback (an alignment).
/// `all` selects whether all end positions with cost <= k are reported or only the rightmost
/// position of each local minimum (sassy's default).
///
/// Returns null if `profile` is invalid.
#[no_mangle]
pub extern "C" fn sassy_capi_new(profile: *const c_char, trace: bool, all: bool) -> *mut Handle {
    if profile.is_null() { return std::ptr::null_mut(); }
    let name = match unsafe { CStr::from_ptr(profile) }.to_str() {
        Ok(s) => s.to_ascii_lowercase(),
        Err(_) => return std::ptr::null_mut(),
    };

    let searcher = match name.as_str() {
        "ascii" => SearcherType::Ascii(if trace { Searcher::<Ascii>::new_fwd() }
                                       else    { Searcher::<Ascii>::new_fwd().without_trace() }),
        "dna"   => SearcherType::Dna(if trace   { Searcher::<Dna>::new_fwd() }
                                     else       { Searcher::<Dna>::new_fwd().without_trace() }),
        "iupac" => SearcherType::Iupac(if trace { Searcher::<Iupac>::new_fwd() }
                                       else     { Searcher::<Iupac>::new_fwd().without_trace() }),
        _ => return std::ptr::null_mut(),
    };

    Box::into_raw(Box::new(Handle { searcher, all }))
}

/// Frees a searcher created by `sassy_capi_new`.
#[no_mangle]
pub extern "C" fn sassy_capi_free(handle: *mut Handle) {
    if handle.is_null() { return; }
    unsafe { drop(Box::from_raw(handle)); }
}

/// Searches `pattern` in `text` with at most `k` edits and returns the number of matches whose end
/// position in `text` is >= `min_end` (the matches ending before `min_end` lie in the overlap with
/// the previous chunk and have already been counted).
///
/// If `out_bytes` is not null, the size in bytes of the match vector returned by sassy is written
/// to it.
#[no_mangle]
pub extern "C" fn sassy_capi_search(
    handle: *mut Handle,
    pattern: *const u8, pattern_len: usize,
    text: *const u8, text_len: usize,
    k: usize, min_end: usize,
    out_bytes: *mut u64,
) -> u64 {
    if handle.is_null() || pattern.is_null() || text.is_null() { return 0; }
    let handle = unsafe { &mut *handle };
    let pattern = unsafe { slice::from_raw_parts(pattern, pattern_len) };
    let text = unsafe { slice::from_raw_parts(text, text_len) };

    macro_rules! run {
        ($s:expr) => { if handle.all { $s.search_all(pattern, &text, k) }
                       else          { $s.search(pattern, &text, k) } };
    }
    let matches = match &mut handle.searcher {
        SearcherType::Ascii(s) => run!(s),
        SearcherType::Dna(s)   => run!(s),
        SearcherType::Iupac(s) => run!(s),
    };

    if !out_bytes.is_null() {
        let bytes = matches.capacity() as u64 * std::mem::size_of::<sassy::Match>() as u64;
        unsafe { *out_bytes = bytes; }
    }

    matches.iter().filter(|m| m.text_end >= min_end).count() as u64
}

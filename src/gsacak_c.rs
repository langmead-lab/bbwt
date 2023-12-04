use std::ptr;

extern "C" {
    fn gsacak(text: *const u8, suffix_array: *mut u32, lcp: *mut i32, da: *mut i32, len: u32) -> i32;

    fn sacak_int(text: *const u32, suffix_array: *mut u32, n: u32, k: u32) -> i32;
}

pub(crate) unsafe fn sacak_u32(text: impl AsRef<[u32]>, suffix_array: &mut [u32]) {
    let text = text.as_ref();
    let k = text.iter().max().unwrap();
    let n = text.len() as u32;

    // we ignore the return value for now
    unsafe{
            sacak_int(text.as_ptr(), suffix_array.as_mut_ptr(), n, *k + 1);

    }
}

pub(crate) unsafe fn gsacak_u8(text: impl AsRef<[u8]>, suffix_array: &mut [u32], lcp: Option<&mut [i32]>) {
    let text = text.as_ref();
    let n = text.len() as u32;
    let lcp = if let Some(arr) = lcp { arr.as_mut_ptr() } else { ptr::null_mut() };

    unsafe {
           gsacak(text.as_ptr(), suffix_array.as_mut_ptr(), lcp, ptr::null_mut(), n);
    }
}

use ducksmiles_analytics::{ds_analytics_free, ds_analytics_run};
use std::ptr;

#[test]
fn ffi_returns_owned_json_and_rejects_shape_overflow_without_dereference() {
    // The buffer length/ownership contract and pre-dereference shape validation matter.
    let op = b"acquisition";
    let data = [2., 1.];
    let args = [2., 0., 0.];
    unsafe {
        let mut len = 0;
        let result = ds_analytics_run(
            op.as_ptr(),
            op.len(),
            data.as_ptr(),
            1,
            2,
            args.as_ptr(),
            3,
            &mut len,
        );
        assert!(!result.is_null());
        let text = std::str::from_utf8(std::slice::from_raw_parts(result, len)).unwrap();
        assert!(text.starts_with("{\"version\":1,"));
        assert!(text.contains("0.398942280"));
        ds_analytics_free(result, len);
        assert!(ds_analytics_run(
            op.as_ptr(),
            op.len(),
            data.as_ptr(),
            usize::MAX,
            2,
            args.as_ptr(),
            3,
            &mut len
        )
        .is_null());
        assert_eq!(len, 0);
        assert!(
            ds_analytics_run(ptr::null(), 0, ptr::null(), 0, 0, ptr::null(), 0, &mut len).is_null()
        );
        ds_analytics_free(ptr::null_mut(), 0);
    }
}

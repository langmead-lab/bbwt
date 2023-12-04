use cc;

fn main() {
    cc::Build::new()
        .flag("-Wno-sign-compare")
        .flag("-Wno-unused-parameter")
        .file("src/gsa/gsacak.c")
        .compile("gsacak");
}

b:
	clear
	cargo build
d:
	clear
	RUSTDOCFLAGS="--html-in-header katex-header.html" cargo doc --no-deps --open
	# https://docs.rs/rustdoc-katex-demo/0.1.5/rustdoc_katex_demo/

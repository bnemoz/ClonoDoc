//! Build script: embed Windows version/publisher metadata into the `.exe`.
//! No-op on non-Windows targets. Failures are non-fatal (metadata is best-effort).

fn main() {
    #[cfg(windows)]
    {
        let mut res = winresource::WindowsResource::new();
        res.set("ProductName", "ClonoDoc");
        res.set("FileDescription", "ClonoDoc — antibody cloning verifier");
        res.set("CompanyName", "Scripps Research");
        res.set("LegalCopyright", "MIT License");
        res.set("OriginalFilename", "clonodoc.exe");
        // Best-effort: if the resource compiler is unavailable, build without metadata.
        let _ = res.compile();
    }
}

# scDLO Hi-C

Single cell DLO Hi-C analysis pipeline.

## Dependency

1. Install [DLO-HiC-Tools](https://github.com/GangCaoLab/DLO-HiC-Tools#installation).
2. Clone this repo, and compile the `extpet` tool.

```bash
# Install Rust
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
$ export PATH=$PATH:~/.cargo/bin  # add cargo path to the PATH env variable
# Compile
$ cd extpet; cargo build --release
```

3. Add `<path_the_repo>/extpet/target/release` and `<path_the_repo>/script` to PATH in your `.bashrc`.


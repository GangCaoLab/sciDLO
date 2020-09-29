# sciDLO Hi-C

Single cell indexed DLO Hi-C analysis pipeline.

## Install

1. Install [DLO-HiC-Tools](https://github.com/GangCaoLab/DLO-HiC-Tools#installation).
2. Install [Nextflow]https://www.nextflow.io/)(version >= 20.07.1.5412)
3. Clone this repo, and compile the `extpet` tool.

```bash
# Install Rust
$ curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
$ export PATH=$PATH:~/.cargo/bin  # add cargo path to the PATH env variable
# Compile
$ cd scTools; cargo build --release
```

4. Add `export PATH=<path_the_repo>/scripts:$PATH` to the end of your `.bashrc`.

## Run pipeline

1. Copy `main.nf` and `nextflow.config` to your working path.
2. Modify `nextflow.config`.
3. Run `nextflow run main.nf`.

## Result

Result file of each steps will stored in the `result` folder in the working path.

| sub folder | description |
| -----------| ------------|
| pairs_cell | pairs of each cell |
| pairs_lib  | pairs of each library |
| pairs_all  | merged pairs |
| dot_hic_cell | .hic file of each cell |
| dot_hic_lib | .hic file of each library |
| dot_hic_all | .hic file of all merged reads |
| tmp | intermediate files of each steps |

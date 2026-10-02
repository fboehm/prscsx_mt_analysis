configfile: "config.yaml"


# Build wildcard combinations
SCENARIOS = expand(
    "{n_eur}_{n_eas}_{rg}_{seed}",
    n_eur=config["n_eur"],
    n_eas=config["n_eas"],
    rg=config["rg"],
    seed=config["seeds"],
)

TRAITS = [0, 1]


rule all:
    input:
        "results/summary.csv",


rule simulate:
    output:
        bim=    "results/sim/{scenario}/data/sim_data.bim",
        refpanel="results/sim/{scenario}/data/snpinfo_mult_1kg_hm3",
    params:
        out_dir="results/sim/{scenario}/data",
        prscsx_mt_dir=config["prscsx_mt_dir"],
        n_snp=config["n_snp"],
        n_causal=config["n_causal"],
        block_size=config["block_size"],
        ld_decay=config["ld_decay"],
        rho_pop=config["rho_pop"],
        h2=config["h2"],
        populations=config["pop"],
        chrom=config["chrom"],
        sim_mode=config.get("sim_mode", "analytical"),
    script:
        "scripts/simulate.py"


rule run_prscsx:
    input:
        bim="results/sim/{scenario}/data/sim_data.bim",
        refpanel="results/sim/{scenario}/data/snpinfo_mult_1kg_hm3",
    output:
        touch("results/prscsx/{scenario}/trait{trait}/output/.done"),
    params:
        data_dir="results/sim/{scenario}/data",
        out_dir="results/prscsx/{scenario}/trait{trait}/output",
        prscsx_mt_dir=config["prscsx_mt_dir"],
        populations=config["pop"],
        chrom=config["chrom"],
        a=config["a"],
        b=config["b"],
        n_iter=config["n_iter"],
        n_burnin=config["n_burnin"],
        thin=config["thin"],
    script:
        "scripts/run_single_trait.py"


rule run_prscsx_mt:
    input:
        bim="results/sim/{scenario}/data/sim_data.bim",
        refpanel="results/sim/{scenario}/data/snpinfo_mult_1kg_hm3",
    output:
        touch("results/prscsx_mt/{scenario}/output/.done"),
    params:
        data_dir="results/sim/{scenario}/data",
        out_dir="results/prscsx_mt/{scenario}/output",
        prscsx_mt_dir=config["prscsx_mt_dir"],
        populations=config["pop"],
        chrom=config["chrom"],
        a=config["a"],
        b=config["b"],
        n_iter=config["n_iter"],
        n_burnin=config["n_burnin"],
        thin=config["thin"],
    script:
        "scripts/run_multi_trait.py"


rule evaluate_prscsx:
    input:
        done="results/prscsx/{scenario}/trait{trait}/output/.done",
        bim="results/sim/{scenario}/data/sim_data.bim",
    output:
        json="results/eval/prscsx/{scenario}/trait{trait}.json",
    params:
        data_dir="results/sim/{scenario}/data",
        out_dir="results/prscsx/{scenario}/trait{trait}/output",
        prscsx_mt_dir=config["prscsx_mt_dir"],
        populations=config["pop"],
        chrom=config["chrom"],
        a=config["a"],
        b=config["b"],
        method="prscsx",
    script:
        "scripts/evaluate.py"


rule evaluate_prscsx_mt:
    input:
        done="results/prscsx_mt/{scenario}/output/.done",
        bim="results/sim/{scenario}/data/sim_data.bim",
    output:
        json="results/eval/prscsx_mt/{scenario}/trait{trait}.json",
    params:
        data_dir="results/sim/{scenario}/data",
        out_dir="results/prscsx_mt/{scenario}/output",
        prscsx_mt_dir=config["prscsx_mt_dir"],
        populations=config["pop"],
        chrom=config["chrom"],
        a=config["a"],
        b=config["b"],
        method="prscsx_mt",
    script:
        "scripts/evaluate.py"


rule aggregate:
    input:
        prscsx=expand(
            "results/eval/prscsx/{scenario}/trait{trait}.json",
            scenario=SCENARIOS,
            trait=TRAITS,
        ),
        prscsx_mt=expand(
            "results/eval/prscsx_mt/{scenario}/trait{trait}.json",
            scenario=SCENARIOS,
            trait=TRAITS,
        ),
    output:
        csv="results/summary.csv",
    script:
        "scripts/aggregate.py"

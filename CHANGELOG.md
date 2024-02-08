# Changelog

## [2.6.0](https://github.com/xsitarcik/amrWatch/compare/v2.5.0...v2.6.0) (2024-02-08)


### Features

* added gtbdtk taxa conversion to ncbi ([3371f0b](https://github.com/xsitarcik/amrWatch/commit/3371f0bf73af2c757d3ee3ecc73cc8b30d735149))

## [2.5.0](https://github.com/xsitarcik/amrWatch/compare/v2.4.0...v2.5.0) (2024-02-08)


### Features

* added db to be set in config for mobsuite as ncbi_plasmids_db_dir ([2bcac0f](https://github.com/xsitarcik/amrWatch/commit/2bcac0f38bc5490d4c2bb88b7043aac75e5fe6ce))
* added SISTR for salmonella ([3b102e1](https://github.com/xsitarcik/amrWatch/commit/3b102e1053f4a1eef5d83641613581e3ba3d06b8))


### Bug Fixes

* fixed rerunning gtdb classification if not empty directory ([46cdde7](https://github.com/xsitarcik/amrWatch/commit/46cdde78a2d154b66793ee8afd67e99b2451bcfe))

## [2.4.0](https://github.com/xsitarcik/amrWatch/compare/v2.3.0...v2.4.0) (2024-02-07)


### Features

* added db download for SCCmec ([9dae50d](https://github.com/xsitarcik/amrWatch/commit/9dae50d2c7a495a1d13acdf5cca2db030322bf25))
* added nextseq trimming mode to cutadapt ([56c7894](https://github.com/xsitarcik/amrWatch/commit/56c7894e6b205db6fdf0f69253fd94247410b221))
* added sccmec ([fd08726](https://github.com/xsitarcik/amrWatch/commit/fd08726106b510e0b1f91f941861ac0d7f79d87c))


### Bug Fixes

* added sccmec for staphyloccocus aureus only ([50e0829](https://github.com/xsitarcik/amrWatch/commit/50e08296653bd914e8d5c39049c7ac3cd8b8bb7b))
* added seqkit stats to count number of contigs ([a6fc872](https://github.com/xsitarcik/amrWatch/commit/a6fc8729879c2334208b356243211bd51e0f324f))
* nextseq trimming is boolean in the config and the cutoff value is taken from r1 3 quality ([e13f2bc](https://github.com/xsitarcik/amrWatch/commit/e13f2bce2b0f6918ad614f39e01e2ce37570563f))

## [2.3.0](https://github.com/xsitarcik/amrWatch/compare/v2.2.1...v2.3.0) (2024-02-06)


### Features

* added mob_suite tools ([9f09fae](https://github.com/xsitarcik/amrWatch/commit/9f09fae937007490cccaad325194840f32b8146f))

## [2.2.1](https://github.com/xsitarcik/amrWatch/compare/v2.2.0...v2.2.1) (2024-02-04)


### Bug Fixes

* added threshold for max assembly length ([43468d8](https://github.com/xsitarcik/amrWatch/commit/43468d891afdb02f1cfafe918eae5a120eb291e8))

## [2.2.0](https://github.com/xsitarcik/amrWatch/compare/v2.1.0...v2.2.0) (2024-01-31)


### Features

* added updateable database for spatyper ([c9af3b2](https://github.com/xsitarcik/amrWatch/commit/c9af3b2d9f32315080e79de53fb66b82e6624e69))

## [2.1.0](https://github.com/xsitarcik/amrWatch/compare/v2.0.2...v2.1.0) (2024-01-29)


### Features

* added abricate parameters for coverage and identity ([114680b](https://github.com/xsitarcik/amrWatch/commit/114680b7c4828725004353d184270547351e3327))
* added check for assembly minimum length ([f9239ad](https://github.com/xsitarcik/amrWatch/commit/f9239ad05ed7aae9bf41994cf8feb738740a24a1))


### Bug Fixes

* forced unicycler to keep only final results ([0a9a135](https://github.com/xsitarcik/amrWatch/commit/0a9a13540538421a2267b326aeac41403b5a068b))
* when restarting gtdbtk classify force remove previous results ([583586d](https://github.com/xsitarcik/amrWatch/commit/583586dfa15fd3bce63890a67b006100f255f06f))

## [2.0.2](https://github.com/xsitarcik/amrWatch/compare/v2.0.1...v2.0.2) (2024-01-23)


### Bug Fixes

* fixed parsing bandage report in case of zero ([2509d07](https://github.com/xsitarcik/amrWatch/commit/2509d07c1a934eeb156b83441e7342e0b121e9e8))
* matching for amrfinder and mlst happens first for genus and species, and then just for genus ignoring entries if they have species ([f837e91](https://github.com/xsitarcik/amrWatch/commit/f837e91f57f8637c9a2d1c44203208ca246801ed))
* mlst and amrfinder receives no organism if not found in map, instead of error ([5b2740e](https://github.com/xsitarcik/amrWatch/commit/5b2740e415a071fd480be4a59ea566924f900a07))

## [2.0.1](https://github.com/xsitarcik/amrWatch/compare/v2.0.0...v2.0.1) (2024-01-22)


### Bug Fixes

* allow to pass also WARN samples ([c7697d6](https://github.com/xsitarcik/amrWatch/commit/c7697d68aa24d6abf322f8bb23618e7792a8b948))
* assembly running mode works correctly per sample ([59a7d23](https://github.com/xsitarcik/amrWatch/commit/59a7d231ce80597287667fae45c4a2ab96cae6a6))

## [2.0.0](https://github.com/xsitarcik/amrWatch/compare/v1.1.0...v2.0.0) (2024-01-21)


### ⚠ BREAKING CHANGES

* reworked amrfinder with input user db

### Features

* added gtdb mapping for amrfinder and mlst ([a5ff773](https://github.com/xsitarcik/amrWatch/commit/a5ff773820b1047ceb0ccbaafad05e1a145c851a))
* added running from assembly ([23d0a34](https://github.com/xsitarcik/amrWatch/commit/23d0a343f40bd77e73a594a56418aba92c3b5a19))
* added taxa dependent outputs, etoki, kleborate and spatyper ([86c59be](https://github.com/xsitarcik/amrWatch/commit/86c59bef973f14ff63c900d5bb8dea97f8267084))
* amrfinder update db automatically ([0108ba1](https://github.com/xsitarcik/amrWatch/commit/0108ba107bf701229a1efed4c1ed61353abcda3c))
* DB for abricate is selectable from config ([14a1c26](https://github.com/xsitarcik/amrWatch/commit/14a1c264c6eba6e2a77b2c156f819823103a1966))
* expanded cutadapt params in config, added mode for assembly input ([2dad26c](https://github.com/xsitarcik/amrWatch/commit/2dad26c2c6456be8101e8c7aae5005f51e8457ca))
* replaced abritamr by amrfinderplus ([1b7d9e8](https://github.com/xsitarcik/amrWatch/commit/1b7d9e847eb95af2dbff954ae2051e1f9dab1110))
* reworked amrfinder with input user db ([d069444](https://github.com/xsitarcik/amrWatch/commit/d0694447fd0ee490fcefdb1dcb119a07d8f50d00))

## [1.1.0](https://github.com/xsitarcik/amrWatch/compare/v1.0.0...v1.1.0) (2024-01-15)


### Features

* added coverage check from qualimap ([03731a2](https://github.com/xsitarcik/amrWatch/commit/03731a2af85228ede82e99688ef354fd0b8a9322))
* added selfcontamination check ([1d65f5f](https://github.com/xsitarcik/amrWatch/commit/1d65f5f4ce4175327a1aa43235b505fe374acff5))
* reworked cutadapt to produce output for multiqc ([919dc17](https://github.com/xsitarcik/amrWatch/commit/919dc17b54f2e84261bdb73a3c10649554532fce))


### Bug Fixes

* added genera check threshold to be set from config ([90bbe18](https://github.com/xsitarcik/amrWatch/commit/90bbe185960acd7be860ed8c4d29b78814b915af))
* fixed cutadapt params for end cutting ([4e7e797](https://github.com/xsitarcik/amrWatch/commit/4e7e79727ac1f2f573e2a0c3b5800ee8f7d16d48))
* multiqc produced for all samples ([9c095a4](https://github.com/xsitarcik/amrWatch/commit/9c095a4ea5083f74b2a4e153e9be3232b12e1726))
* multiqc trims kreport2 from sample ([192719d](https://github.com/xsitarcik/amrWatch/commit/192719db56bc20314d7153823a451915e7f58b21))

## 1.0.0 (2024-01-14)


### Features

* first version of amrwatch ([0ce9349](https://github.com/xsitarcik/amrWatch/commit/0ce9349d055e72eae97999da1a68614399568af9))

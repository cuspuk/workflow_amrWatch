# Changelog

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

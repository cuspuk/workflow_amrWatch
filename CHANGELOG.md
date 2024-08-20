# Changelog

## [2.21.4](https://github.com/cuspuk/workflow_amrWatch/compare/v2.21.3...v2.21.4) (2024-08-20)


### Bug Fixes

* added assembly construction check to summary ([1cc655a](https://github.com/cuspuk/workflow_amrWatch/commit/1cc655af6954b883b85a411fa41ae11b410db0c5))
* empty parsed taxa fixed for resfinder ([7e64085](https://github.com/cuspuk/workflow_amrWatch/commit/7e64085fd41b01262aa0fb68083501147ee53693))

## [2.21.3](https://github.com/cuspuk/workflow_amrWatch/compare/v2.21.2...v2.21.3) (2024-08-09)


### Bug Fixes

* ncbi_gtdb_mappings ([d2463c9](https://github.com/cuspuk/workflow_amrWatch/commit/d2463c9e1b73f4c9d3ec0f3041d4652c410628de))

## [2.21.2](https://github.com/cuspuk/workflow_amrWatch/compare/v2.21.1...v2.21.2) (2024-08-05)


### Bug Fixes

* reverted gzipping fasta ([c29b643](https://github.com/cuspuk/workflow_amrWatch/commit/c29b6437743050b91388006ddac4e54efff442da))

## [2.21.1](https://github.com/cuspuk/workflow_amrWatch/compare/v2.21.0...v2.21.1) (2024-08-01)


### Bug Fixes

* added gzipping of fasta ([fdd174e](https://github.com/cuspuk/workflow_amrWatch/commit/fdd174e22ca6211695a2b9f5341be8b6cc273d2f))

## [2.21.0](https://github.com/cuspuk/workflow_amrWatch/compare/v2.20.3...v2.21.0) (2024-07-27)


### Features

* expanded foreign contamination check to first merge similar groups ([9ae6fae](https://github.com/cuspuk/workflow_amrWatch/commit/9ae6fae5830d2b66d12f424546496c8c8aec40ff))

## [2.20.3](https://github.com/cuspuk/workflow_amrWatch/compare/v2.20.2...v2.20.3) (2024-07-26)


### Bug Fixes

* added custom mapping to ncbi ID to the config ([27b034d](https://github.com/cuspuk/workflow_amrWatch/commit/27b034d5dd062f596fb9d7c2d516f52fa40fc9b6))
* reworked ncbi id conversion to produce null if ambiguous ([69bcd5b](https://github.com/cuspuk/workflow_amrWatch/commit/69bcd5b1e238477a21711e74c606af49340d3b06))

## [2.20.2](https://github.com/cuspuk/workflow_amrWatch/compare/v2.20.1...v2.20.2) (2024-07-25)


### Bug Fixes

* number of basepairs in check now calculated from seqkit stats ([82b8377](https://github.com/cuspuk/workflow_amrWatch/commit/82b837709718fd9cc87ebc306728eeaeafb82a12))

## [2.20.1](https://github.com/cuspuk/workflow_amrWatch/compare/v2.20.0...v2.20.1) (2024-06-28)


### Bug Fixes

* ncbi taxid parsing ([5527eb1](https://github.com/cuspuk/workflow_amrWatch/commit/5527eb104dd52fe02fecd85b611a7225a0e949da))

## [2.20.0](https://github.com/cuspuk/workflow_amrWatch/compare/v2.19.0...v2.20.0) (2024-06-28)


### Features

* added bam stats using samtools and picard insert ([4debace](https://github.com/cuspuk/workflow_amrWatch/commit/4debace319fe66ae3bfd03fc4058852ddd1de026))
* added new mapping for new mlst db version ([71168fc](https://github.com/cuspuk/workflow_amrWatch/commit/71168fc3ac45faabb9f0e7ae612c7d0e08fd7967))
* assembly statistics using seqkit all ([a37ba5e](https://github.com/cuspuk/workflow_amrWatch/commit/a37ba5eafae85e8190c55e1e685ee27be4c1f0ed))
* bumped modules and gtdbtk ([7680fc9](https://github.com/cuspuk/workflow_amrWatch/commit/7680fc918aae99be5fb39a374ac1c4d4ca6e02aa))


### Bug Fixes

* final bam and bai kept ([d8c6b77](https://github.com/cuspuk/workflow_amrWatch/commit/d8c6b7784d02c1b7e437d3881db4d9c7f2199797))

## [2.19.0](https://github.com/cuspuk/workflow_amrWatch/compare/v2.18.2...v2.19.0) (2024-06-20)


### Features

* added amrfinder only output summary ([98a7e22](https://github.com/cuspuk/workflow_amrWatch/commit/98a7e22d75b2609938010f2a880ade1ce570c1b8))

## [2.18.2](https://github.com/cuspuk/workflow_amrWatch/compare/v2.18.1...v2.18.2) (2024-05-22)


### Bug Fixes

* added escherichia_o_antigen, escherichia_h_antigen and pneumo_capsular_type into summary ([76a0cde](https://github.com/cuspuk/workflow_amrWatch/commit/76a0cdeeb31814a63476e0a56960953bfe20e2f7))
* added ncbi_taxonomy_id into summary ([376d7d5](https://github.com/cuspuk/workflow_amrWatch/commit/376d7d543bdd0375b04a77cdc1c5bad7e0d860c3))
* added serotype column filled by either seqsero, ebeis or pneumokity ([5e6063f](https://github.com/cuspuk/workflow_amrWatch/commit/5e6063f2ad44c0e5f9d85b6e80bf7ce7e24ce6f8))
* ecoli and pneumokity serotype columns are added after salmonella columns ([96b2efc](https://github.com/cuspuk/workflow_amrWatch/commit/96b2efc4a8df6913f394f62fec399a715737442c))
* handled json string in csv ([cbbbff1](https://github.com/cuspuk/workflow_amrWatch/commit/cbbbff128728d1025328ea4f48af5f6cbe2fbb2a))
* pneumokitty outputs fixed to produce a new dir in amr_detect folder ([fb28967](https://github.com/cuspuk/workflow_amrWatch/commit/fb28967cfb0f31e8f0a57d65bad098ce057d6848))
* pneumokity download handles already existing dir ([5dbd3fd](https://github.com/cuspuk/workflow_amrWatch/commit/5dbd3fdf0d0796fd65736f39e232348c17f558ca))


### Performance Improvements

* bumped pneumokity down to python3.8 ([4d6cc31](https://github.com/cuspuk/workflow_amrWatch/commit/4d6cc31ce42609ebb82439c6dad22fc424bc6408))

## [2.18.1](https://github.com/cuspuk/workflow_amrWatch/compare/v2.18.0...v2.18.1) (2024-05-15)


### Bug Fixes

* changes after transfer repo ([8a85f79](https://github.com/cuspuk/workflow_amrWatch/commit/8a85f79469baf35fb40d9111b4f19656f5e0acfe))

## [2.18.0](https://github.com/xsitarcik/amrWatch/compare/v2.17.0...v2.18.0) (2024-05-13)


### Features

* added pneumokity for Pneumococcus ([ca91929](https://github.com/xsitarcik/amrWatch/commit/ca91929ac76b64ebe0f02c6f66bf4bd455d5952f))


### Bug Fixes

* pneumokity for Streptococcus pneumoniae ([f16f158](https://github.com/xsitarcik/amrWatch/commit/f16f1585a836bdbbf6bc2cb5ee6be243bdd26ef0))

## [2.17.0](https://github.com/xsitarcik/amrWatch/compare/v2.16.2...v2.17.0) (2024-04-19)


### Features

* added insilico PCR ([d31dbe7](https://github.com/xsitarcik/amrWatch/commit/d31dbe79d306e8a1eb347e41b4e2c5763664cc98))


### Bug Fixes

* added max_mismatch parameter for insilico ([eb46309](https://github.com/xsitarcik/amrWatch/commit/eb4630981718530c26af7646cfd983f5a210d342))
* downgrade seqkit for amplicon due to bug ([dde81fb](https://github.com/xsitarcik/amrWatch/commit/dde81fb2c406283e34344da230f22a53395eff74))

## [2.16.2](https://github.com/xsitarcik/amrWatch/compare/v2.16.1...v2.16.2) (2024-04-11)


### Bug Fixes

* failed preassembly check fixed to not request assembly ([4f5d57d](https://github.com/xsitarcik/amrWatch/commit/4f5d57d3dc5f433589293a8e07a57675507d8c36))

## [2.16.1](https://github.com/xsitarcik/amrWatch/compare/v2.16.0...v2.16.1) (2024-03-28)


### Bug Fixes

* added flag for controlling hamronization running ([62749aa](https://github.com/xsitarcik/amrWatch/commit/62749aad07f181a3e5367b6c7e7b6d88f803a95c))
* hamronization flag excludes even running per sample ([edf4e40](https://github.com/xsitarcik/amrWatch/commit/edf4e40f39e4b498dfc4eeb7a4a5ce7d715cb453))

## [2.16.0](https://github.com/xsitarcik/amrWatch/compare/v2.15.2...v2.16.0) (2024-03-27)


### Features

* added clonal complex ([232aab0](https://github.com/xsitarcik/amrWatch/commit/232aab07c332c080f090b0e1dea04fc0205ee980))


### Bug Fixes

* clonal complex set to dash if found but empty ([5d35833](https://github.com/xsitarcik/amrWatch/commit/5d35833c73347b48079ccb163c91757b8cc1dc67))
* clonal complex set to dash if not found in profile ([1a82d49](https://github.com/xsitarcik/amrWatch/commit/1a82d49b1024184bdeb5647121fc8bf76294b5f3))
* schema requires ccc ([9b2da8d](https://github.com/xsitarcik/amrWatch/commit/9b2da8d1243b19943e2f04310b8dd7303c90637b))

## [2.15.2](https://github.com/xsitarcik/amrWatch/compare/v2.15.1...v2.15.2) (2024-03-27)


### Bug Fixes

* fixed spatyper in summary, added sccmec output into summary ([ad5a039](https://github.com/xsitarcik/amrWatch/commit/ad5a039520be699f35f42f6f05a1de757e43c3c5))
* sccmec added merged column with only true types ([9618725](https://github.com/xsitarcik/amrWatch/commit/96187255488755635d0f04dca236be7e0dc0702d))
* sccmec also in merged summary ([e1ef7af](https://github.com/xsitarcik/amrWatch/commit/e1ef7af711fb535d4fab4dbe3181ef783077fe15))

## [2.15.1](https://github.com/xsitarcik/amrWatch/compare/v2.15.0...v2.15.1) (2024-03-26)


### Bug Fixes

* fixed logging for mlst skipping ([4e06cb8](https://github.com/xsitarcik/amrWatch/commit/4e06cb8e3a20b0afbe44b79c16cd8ee0622a0b43))
* reworked main workflow logs to use main snakemake logger ([23e2663](https://github.com/xsitarcik/amrWatch/commit/23e2663b0a16aa84ca51f2d771057ddb93c2f365))


### Performance Improvements

* fixed some logging statements ([4c9ee1c](https://github.com/xsitarcik/amrWatch/commit/4c9ee1c9c1ee14001f790dcc6784510ae4049358))

## [2.15.0](https://github.com/xsitarcik/amrWatch/compare/v2.14.2...v2.15.0) (2024-03-25)


### Features

* added assembly__unicycler-&gt;use_long_if_relevant element in config to ignore long reads even if provided ([929cb78](https://github.com/xsitarcik/amrWatch/commit/929cb782733ab4b2eae87116b11fdc5daa5f6621))
* added long column in samples.csv, which automatically runs unicycler with long reads ([cf429c9](https://github.com/xsitarcik/amrWatch/commit/cf429c942cdf1307d71c683175e593577d2ccf3a))


### Bug Fixes

* abricate values are sorted in summary by human sorting ([9227ffb](https://github.com/xsitarcik/amrWatch/commit/9227ffb018b5e5f68bb91a4d98a6bd1b0927da09))
* abricate virulence hits sort is now case insensitive ([82b6b25](https://github.com/xsitarcik/amrWatch/commit/82b6b25ee4c35d614598a49d6e7f5d314dc03dc9))
* assembly fasta headers cleaned, added prefix "_sample_contig_", replaced spaces with double underscores, and equal sign with single underscore ([c03a88f](https://github.com/xsitarcik/amrWatch/commit/c03a88f867172f6f56fd83a07b8b55e817e99f82))
* assembly with cleaned headers is now used for subsequent steps instead of original fasta assembly ([d5b0847](https://github.com/xsitarcik/amrWatch/commit/d5b0847e8307657776c6c2680d33cb52cd38bbd4))
* each qc_check is represented in summary by three columns {check}__value, {check}__result, {check}__comment ([aeca81b](https://github.com/xsitarcik/amrWatch/commit/aeca81b8dc3b43f3ff3af06343a73f3567515fad))
* log paths changed for resfinder db ([1344a3a](https://github.com/xsitarcik/amrWatch/commit/1344a3aeb4e616bbd592cad1006851718b960170))
* mlst is requested only for taxa with matched scheme ([4da5352](https://github.com/xsitarcik/amrWatch/commit/4da535231be4c1a006fa587665b6bc7d64cfc1ba))
* removed leading underscore before sample name in fasta header ([0ddb2bb](https://github.com/xsitarcik/amrWatch/commit/0ddb2bb99759c5fce0cb05e2c14231a4af4537b1))
* reordered qc checks to result, value and comment, for each column ([2813a10](https://github.com/xsitarcik/amrWatch/commit/2813a1079a9e3a877d18d21c25beba092c415b43))

## [2.14.2](https://github.com/xsitarcik/amrWatch/compare/v2.14.1...v2.14.2) (2024-03-17)


### Bug Fixes

* added kma_index needed for resfinder ([1851d8e](https://github.com/xsitarcik/amrWatch/commit/1851d8e8a8426fb996f2170fd9b2747ce0626bc6))
* hamronize skips resfinder if fastq ([a755088](https://github.com/xsitarcik/amrWatch/commit/a7550889c3766daf2a98a3a2238540e32a9e4ba4))
* skip resfinder even for per sample summary hamronization ([ba660c3](https://github.com/xsitarcik/amrWatch/commit/ba660c3fee56f1e5489f6ac43d5f617248a3fd63))

## [2.14.1](https://github.com/xsitarcik/amrWatch/compare/v2.14.0...v2.14.1) (2024-03-17)


### Bug Fixes

* added assembly_length into summary ([88e3707](https://github.com/xsitarcik/amrWatch/commit/88e370778e0a22b58e2a98766c1d41540c089dfc))
* seqkit parsed into summary ([6b2c835](https://github.com/xsitarcik/amrWatch/commit/6b2c83593343c398af2e8d278d416ebe85ea0365))

## [2.14.0](https://github.com/xsitarcik/amrWatch/compare/v2.13.0...v2.14.0) (2024-03-16)


### Features

* added merge summary for all samples ([119bddd](https://github.com/xsitarcik/amrWatch/commit/119bddd29b962977a3143513984e087c5ed6aa6c))
* added summarization of results per sample ([98de13d](https://github.com/xsitarcik/amrWatch/commit/98de13d5a6c9a1fcfcc46af1b5bae82dcefc7049))


### Bug Fixes

* added reorder for amrfinder columns in summary ([3ebe858](https://github.com/xsitarcik/amrWatch/commit/3ebe858172767f7db3d0ac2b7ef31b86e57c4758))
* check resfinder input validity before workflow start ([4b6a7d8](https://github.com/xsitarcik/amrWatch/commit/4b6a7d828cb9e9c4cb69a55f958c2ff5e772590c))
* human check contamination fixed ([d29c2a0](https://github.com/xsitarcik/amrWatch/commit/d29c2a07c264092c1db8db8225f6bb9febccb3b8))
* merged summary for amrfinder with class having __ as separator ([61f16bd](https://github.com/xsitarcik/amrWatch/commit/61f16bdd2aee226baf7fe692bb8229795248f9dc))
* summary skips seqkit if input are assemblies ([f5e0e28](https://github.com/xsitarcik/amrWatch/commit/f5e0e28e78ce745d0431810f228ef9bb5b4f763f))

## [2.13.0](https://github.com/xsitarcik/amrWatch/compare/v2.12.0...v2.13.0) (2024-03-15)


### Features

* added resfinder option to be run from reads or assembly ([6f04c9d](https://github.com/xsitarcik/amrWatch/commit/6f04c9d20a3b6daa29c29065901344de5f063bdb))

## [2.12.0](https://github.com/xsitarcik/amrWatch/compare/v2.11.0...v2.12.0) (2024-03-15)


### Features

* added check for human contamination ([9d5f819](https://github.com/xsitarcik/amrWatch/commit/9d5f819bc3e975f77ee09c6ab63d17e5e902639d))

## [2.11.0](https://github.com/xsitarcik/amrWatch/compare/v2.10.2...v2.11.0) (2024-03-15)


### Features

* added hamronization ([087efdb](https://github.com/xsitarcik/amrWatch/commit/087efdbb37bdda5624375764c7b2693e5f3150d0))
* added resfinder and pointfinder to hamronization, added hamronization summary on sample basis ([811d5f4](https://github.com/xsitarcik/amrWatch/commit/811d5f4caba217ec8946c11cd4cbfe207d240bb5))
* added version storing for rgi,amrfinder,abricate ([610e636](https://github.com/xsitarcik/amrWatch/commit/610e636c55f7837b3d83a758a07f98d859ad0746))


### Bug Fixes

* added resfinder and pointfinder correct index for hamronization ([539119d](https://github.com/xsitarcik/amrWatch/commit/539119d9a424bd1db84d607449ea1821ef8c9e4b))
* dont request hamronize summary when no samples passed qc ([22d91a1](https://github.com/xsitarcik/amrWatch/commit/22d91a1e24f463393106012a6200d9a33ce68bf0))
* hamronize amrfinder replaced by amrfinderplus ([2625b1a](https://github.com/xsitarcik/amrWatch/commit/2625b1af354999f11e679d64107bf040dfb32ca9))
* hamronize fixed for abricate ([d4e1e59](https://github.com/xsitarcik/amrWatch/commit/d4e1e59538aaa3c71a92d4f835c492b10eb4530a))
* handle case when pointfinder has no output ([5346bc7](https://github.com/xsitarcik/amrWatch/commit/5346bc7701661830570aeba3ddb02889cfde90c9))
* logging for version getters ([bda7f67](https://github.com/xsitarcik/amrWatch/commit/bda7f6786f5960e1e83d0c5f06b2e58d8e53b007))

## [2.10.2](https://github.com/xsitarcik/amrWatch/compare/v2.10.1...v2.10.2) (2024-03-13)


### Performance Improvements

* bumped wrapper versions ([eb86fee](https://github.com/xsitarcik/amrWatch/commit/eb86fee1c7e08fa9ca622ab299fcdc73b200b6b4))

## [2.10.1](https://github.com/xsitarcik/amrWatch/compare/v2.10.0...v2.10.1) (2024-03-12)


### Bug Fixes

* assembly fasta are checked for existence if using from pepfile ([bdf29d8](https://github.com/xsitarcik/amrWatch/commit/bdf29d859c97074fa63b74bf669db656e66c2ee5))

## [2.10.0](https://github.com/xsitarcik/amrWatch/compare/v2.9.0...v2.10.0) (2024-02-24)


### Features

* added iscrispol for salmonella enterica ([9903705](https://github.com/xsitarcik/amrWatch/commit/9903705837c27e46a618a1fd494d47a109a2711d))


### Bug Fixes

* renamed min_bp to trimmed_fastq_length to increase readability ([367b48d](https://github.com/xsitarcik/amrWatch/commit/367b48d36fd778b557009409b85a00aefa7519d8))
* statements about missing taxa in dictionaries are moved to stderr ([e6892ee](https://github.com/xsitarcik/amrWatch/commit/e6892eee99fdc576392ae6ad858a1b735e50d37c))

## [2.9.0](https://github.com/xsitarcik/amrWatch/compare/v2.8.2...v2.9.0) (2024-02-20)


### Features

* added check for bp before assembly, bracken added to be other condition for assembly ([25763b3](https://github.com/xsitarcik/amrWatch/commit/25763b35969369ab360def7f402b0ccdee9a6972))
* reworked checks to produce tsv format to increase parsability ([c6ab8fc](https://github.com/xsitarcik/amrWatch/commit/c6ab8fc6811e7ab3c193c79cdced23558d18d413))


### Bug Fixes

* added check for adapters fasta as it cant be directly checked from snakemake rule ([c581821](https://github.com/xsitarcik/amrWatch/commit/c58182195084e73e8765c43c713f537615f4d6f3))
* added header for qc_summary tsv ([8475964](https://github.com/xsitarcik/amrWatch/commit/8475964702cd1d0575e13c7fa5fc5f89e6767071))
* check tsvs are removed as they are summarized in other file ([6381e45](https://github.com/xsitarcik/amrWatch/commit/6381e45d00b741454d60a6d9074bc79bd455add8))
* fix download for gtdb metadata ([8b21a46](https://github.com/xsitarcik/amrWatch/commit/8b21a4600804ab8882fcd8f044c2bac2a141c78d))
* marked bwa_index in selfcontamination check as temporary ([21930e8](https://github.com/xsitarcik/amrWatch/commit/21930e89ca5b76962ace0646453f1761cb8b654c))
* multiqc now works with bracken instead of kraken ([360e3b1](https://github.com/xsitarcik/amrWatch/commit/360e3b19c381290818f3614dfdc713fbdd40a699))
* request qc outputs even when one sample is provided ([2d8cc30](https://github.com/xsitarcik/amrWatch/commit/2d8cc3049de76c6541f8400e8b93f6190215cf8c))

## [2.8.2](https://github.com/xsitarcik/amrWatch/compare/v2.8.1...v2.8.2) (2024-02-16)


### Bug Fixes

* requested fastqc and bracken outputs when not requesting multiqc ([871ee74](https://github.com/xsitarcik/amrWatch/commit/871ee748ee057acf3fdeeaac1ad617c67eb80f02))

## [2.8.1](https://github.com/xsitarcik/amrWatch/compare/v2.8.0...v2.8.1) (2024-02-11)


### Bug Fixes

* cleanup localDB for RGI using onsuccess instead of temp flag as snakemake cannot handle checkpoints in this case ([14c8dce](https://github.com/xsitarcik/amrWatch/commit/14c8dce6d2ffcdbb328d325df5c198ff2f17c210))

## [2.8.0](https://github.com/xsitarcik/amrWatch/compare/v2.7.0...v2.8.0) (2024-02-10)


### Features

* added check level for self contamination to control whether to warn or to fail ([3a5d26f](https://github.com/xsitarcik/amrWatch/commit/3a5d26f6046a57e7708aca830bed2aa4db807985))
* added resfinder ([d01cc0a](https://github.com/xsitarcik/amrWatch/commit/d01cc0ab9557283796ef5e87751a593c95ffe47f))
* added seqsero2 for salmonella ([99f4075](https://github.com/xsitarcik/amrWatch/commit/99f40754d15043d06ffa090f5597c08371c059d4))
* isolated mlst db and exposed path to be set in config ([53fe1c4](https://github.com/xsitarcik/amrWatch/commit/53fe1c44893a679be87c81df03ed70e8105ee160))


### Bug Fixes

* sistr output produced with tab extension ([bd71abe](https://github.com/xsitarcik/amrWatch/commit/bd71abe81b41f5f05a02376b01877da2ebb8104a))


### Performance Improvements

* bumped sistr_cmd ([5b35a2e](https://github.com/xsitarcik/amrWatch/commit/5b35a2ea41a56990857c0d5c8deaf9c4a00832fa))

## [2.7.0](https://github.com/xsitarcik/amrWatch/compare/v2.6.0...v2.7.0) (2024-02-09)


### Features

* added RGI database and caller ([6e4c3a6](https://github.com/xsitarcik/amrWatch/commit/6e4c3a6ede9bee8c15fc9c453da8aa4ebb2edbdb))

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

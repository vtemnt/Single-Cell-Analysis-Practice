# [POSTECH Internship] Single-cell RNA-seq Analysis Pipeline
# Project: Single-Cell-Analysis-Practice
# ---------------------------------------------------------

# 1. 라이브러리 로드 (Library Load)
library(Seurat)

# 2. 데이터 불러오기 (Data Loading)
# 현재 위치(02_Projects/...)에서 상위로 두 번 나가서 01_RawData에 접근합니다.
pbmc.data <- Read10X_h5("../../01_RawData/pbmc_1k_v3_filtered_feature_bc_matrix.h5")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "POSTECH_Practice")

# 3. 품질 관리 (QC: Quality Control)
# 미토콘드리아 유전자 비율 계산 (Mitochondrial gene percentage)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 유전자 수(200~2500) 및 미토콘드리아 비율(5% 미만) 기준으로 필터링
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 4. 데이터 정규화 (Normalization)
# 세포별 RNA 총량을 동일하게 맞추고 로그 변환 수행
pbmc <- NormalizeData(pbmc)

# 5. 고변동 유전자 선정 (Variable Feature Selection)
# 세포 간 발현 차이가 큰 상위 2000개 유전자를 선정
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 6. 결과 확인 및 완료 메시지
print("--- [Success] Seurat Analysis Pipeline Completed up to Variable Features ---")
print(pbmc)

# 고변동 유전자 상위 10개 출력
top10 <- head(VariableFeatures(pbmc), 10)
print("Top 10 Variable Genes:")
print(top10)

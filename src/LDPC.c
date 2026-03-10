
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================= RNG 與 normal() ==================
// 常數（Numerical Recipes 的 ran1）
#define IA   16807L
#define IM   2147483647L
#define AM   (1.0/IM)
#define IQ   127773L
#define IR   2836L
#define NTAB 32
#define NDIV (1 + (IM-1)/NTAB)
#define EPS  1.2e-7
#define RNMX (1.0 - EPS)

// 建立 CSV 檔頭（若不存在）
static void ensure_csv_header(const char *path) {
    FILE *f = fopen(path, "r");
    if (f) { fclose(f); return; }               // 已有就不動
    f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "Algo,SNR_dB,DecodedBits,BitErrors,BER,BlocksTried,BLER,Sec\n");
    fclose(f);
}

// ran1：回傳 (0,1) 開區間的均勻亂數。idum<=0 會初始化。
static double ran1(long *idum) {
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0) *idum += IM;
    j = (int)(iy / NDIV);
    iy = iv[j];
    iv[j] = *idum;
    temp = AM * iy;
    return (temp > RNMX) ? RNMX : temp;
}

// normal： 一次輸出兩個 N(0, σ^2)
static void normal(double *n1, double *n2, double sigma, long *idum) {
    double x1, x2, s;
    do {
        x1 = ran1(idum);
        x2 = ran1(idum);
        x1 = 2.0 * x1 - 1.0;
        x2 = 2.0 * x2 - 1.0;
        s  = x1 * x1 + x2 * x2;
    } while (s >= 1.0 || s == 0.0);
    double factor = sigma * sqrt(-2.0 * log(s) / s);
    *n1 = x1 * factor;
    *n2 = x2 * factor;
}

// ============================= 參數與檔名 =====================================
enum { N = 1023, K = 781, DEG = 32 };

static const double LLR_CLAMP = 50.0;

static const char *FN_SIM = "Sim.txt";
static const char *FN_G   = "ldpc_G_1023.txt";
static const char *FN_H   = "ldpc_H_1023.txt";

// ----  block 的早停門檻 ----
static const unsigned long long STOP_BLOCK_ERRORS = 100ULL;  // 錯 100 個 block 即停
static const unsigned long long STOP_BIT_ERRORS   = 0ULL;    // 若也要「bit 錯 100」就改為 100

// ============================= 小工具 =========================================
static inline double clamp_val(double x, double lo, double hi) {
    return (x < lo) ? lo : ((x > hi) ? hi : x);
}

static inline double clamp_llr(double L) {
    if (L < -LLR_CLAMP) return -LLR_CLAMP;
    if (L >  LLR_CLAMP) return  LLR_CLAMP;
    return L;
}

// 讀取一行僅保留 '0'/'1'（給 G）
static int read_line_01(FILE *fp, char *buf, size_t buf_sz, int *len_out) {
    char line[8192];
    if (!fgets(line, sizeof(line), fp)) return 0;
    int L = 0;
    for (char *p = line; *p; ++p) {
        if (*p == '0' || *p == '1') {
            if (L + 1 < (int)buf_sz) buf[L++] = *p;
        }
    }
    if (len_out) *len_out = L;
    return 1;
}

// ============================= LFSR (m-seq) ==================================
static void make_information_bits(uint8_t *u_all, int64_t blocks) {
    uint8_t s[6] = {1,0,0,0,0,0};
    int64_t total = blocks * (int64_t)K;
    for (int64_t i = 0; i < total; ++i) {
        uint8_t next = s[1] ^ s[0];
        s[0] = s[1];
        s[1] = s[2];
        s[2] = s[3];
        s[3] = s[4];
        s[4] = s[5];
        s[5] = next;
        u_all[i] = s[5] & 1u;
    }
}

// ============================= 讀取 G（K×N） ==================================
static int load_G(uint8_t ***G, int **info_cols) {
    FILE *fp = fopen(FN_G, "r");
    if (!fp) {
        fprintf(stderr, "[ERR] open %s\n", FN_G);
        return 0;
    }
    uint8_t **Gr = (uint8_t**)malloc(sizeof(uint8_t*) * K);
    if (!Gr) { fclose(fp); return 0; }
    for (int i = 0; i < K; ++i) {
        Gr[i] = (uint8_t*)malloc(sizeof(uint8_t) * N);
        if (!Gr[i]) { fclose(fp); return 0; }
    }
    char buf[N + 8];
    int L = 0;
    for (int i = 0; i < K; ++i) {
        if (!read_line_01(fp, buf, sizeof(buf), &L)) {
            fprintf(stderr, "[ERR] G rows不足\n");
            fclose(fp);
            return 0;
        }
        if (L != N) {
            fprintf(stderr, "[ERR] G第%d列長度=%d != %d\n", i, L, N);
            fclose(fp);
            return 0;
        }
        for (int j = 0; j < N; ++j) {
            Gr[i][j] = (buf[j] == '1') ? 1u : 0u;
        }
    }
    fclose(fp);

    // 嘗試找出 systematic 的資訊位欄位（one-hot）
    int *icol = (int*)malloc(sizeof(int) * K);
    if (!icol) return 0;
    for (int i = 0; i < K; ++i) icol[i] = -1;
    int *col_owner = (int*)malloc(sizeof(int) * N);
    if (!col_owner) return 0;
    for (int j = 0; j < N; ++j) col_owner[j] = -1;

    int found = 0;
    for (int j = 0; j < N && found < K; ++j) {
        int cnt = 0, rpos = -1;
        for (int i = 0; i < K; ++i) if (Gr[i][j]) { cnt++; rpos = i; }
        if (cnt == 1 && col_owner[j] == -1 && icol[rpos] == -1) {
            col_owner[j] = rpos;
            icol[rpos]   = j;
            found++;
        }
    }
    if (found != K) {
        fprintf(stderr, "[WRN] 未完整偵測 one-hot；假設資訊位在前 K 欄\n");
        for (int i = 0; i < K; ++i) icol[i] = i;
    }
    free(col_owner);

    *G = Gr;
    *info_cols = icol;
    return 1;
}

// 編碼：x = u G (GF(2))
static void encode_block(uint8_t **G, const uint8_t *u, uint8_t *x) {
    for (int j = 0; j < N; ++j) {
        int acc = 0;
        for (int i = 0; i < K; ++i) {
            acc ^= ((u[i] & 1) & (G[i][j] & 1));
        }
        x[j] = (uint8_t)(acc & 1);
    }
}

// ============================= 讀取 H（鄰接表） ===============================
static int load_H(int **check_to_vars, int **var_to_checks, int **var_pos_in_check) {
    FILE *fp = fopen(FN_H, "r");
    if (!fp) {
        fprintf(stderr, "[ERR] open %s\n", FN_H);
        return 0;
    }
    int total = N * DEG;
    int *CTV  = (int*)malloc(sizeof(int) * total);
    int *VTC  = (int*)malloc(sizeof(int) * total);
    int *VPIC = (int*)malloc(sizeof(int) * total);
    if (!CTV || !VTC || !VPIC) { fclose(fp); return 0; }

    // 前半：check→vars（1-based 轉 0-based）
    for (int m = 0; m < N; ++m) {
        for (int d = 0; d < DEG; ++d) {
            int x = 0;
            if (fscanf(fp, "%d", &x) != 1) {
                fprintf(stderr, "[ERR] H前半讀取不足\n");
                fclose(fp);
                return 0;
            }
            if (x < 1 || x > N) {
                fprintf(stderr, "[ERR] H index超界 %d\n", x);
                fclose(fp);
                return 0;
            }
            CTV[m * DEG + d] = x - 1;
        }
    }
    // 後半：var→checks（1-based 轉 0-based）
    for (int v = 0; v < N; ++v) {
        for (int d = 0; d < DEG; ++d) {
            int x = 0;
            if (fscanf(fp, "%d", &x) != 1) {
                fprintf(stderr, "[ERR] H後半讀取不足\n");
                fclose(fp);
                return 0;
            }
            if (x < 1 || x > N) {
                fprintf(stderr, "[ERR] H index超界 %d\n", x);
                fclose(fp);
                return 0;
            }
            VTC[v * DEG + d] = x - 1;
        }
    }
    fclose(fp);

    // 交叉確認：求 v 在 m 裡的 t 位置（VPIC）
    for (int v = 0; v < N; ++v) {
        for (int s = 0; s < DEG; ++s) {
            int m = VTC[v * DEG + s];
            int t_found = -1;
            for (int t = 0; t < DEG; ++t) {
                if (CTV[m * DEG + t] == v) { t_found = t; break; }
            }
            if (t_found < 0) {
                fprintf(stderr, "[ERR] 不一致：check %d 找不到 var %d\n", m, v);
                return 0;
            }
            VPIC[v * DEG + s] = t_found;
        }
    }

    *check_to_vars = CTV;
    *var_to_checks = VTC;
    *var_pos_in_check = VPIC;
    return 1;
}

// ============================= BPSK 與 σ² ====================================
static inline double bit_to_bpsk(uint8_t b) {
    return (b == 0) ? +1.0 : -1.0;
}

static inline double ebn0_db_to_sigma2(double ebn0_db) {
    double R = (double)K / (double)N;
    double ebn0 = pow(10.0, ebn0_db / 10.0);
    return 1.0 / (2.0 * R * ebn0);
}

// ============================= 解碼資料結構 ===================================
typedef struct {
    int *check_to_vars;    // N*DEG : m*DEG+t -> v
    int *var_to_checks;    // N*DEG : v*DEG+s -> m
    int *var_pos_in_check; // N*DEG : v*DEG+s -> t' (v 在 m 的位置)
    double *L;             // 通道 LLR [N]
    double *q;             // var->check [N*DEG]
    double *r;             // check->var [N*DEG]
} Decoder;

static int decoder_init(Decoder *D, int *CTV, int *VTC, int *VPIC) {
    int total = N * DEG;
    D->check_to_vars = CTV;
    D->var_to_checks = VTC;
    D->var_pos_in_check = VPIC;
    D->L = (double*)malloc(sizeof(double) * N);
    D->q = (double*)malloc(sizeof(double) * total);
    D->r = (double*)malloc(sizeof(double) * total);
    if (!D->L || !D->q || !D->r) {
        fprintf(stderr, "[ERR] OOM decoder\n");
        return 0;
    }
    memset(D->L, 0, sizeof(double) * N);
    memset(D->q, 0, sizeof(double) * total);
    memset(D->r, 0, sizeof(double) * total);
    return 1;
}

static inline int E(int m, int t) {
    return m * DEG + t;
}

static void initialize_llr_and_edges(Decoder *D, const double *y, double sigma2) {
    for (int v = 0; v < N; ++v) {
        double Lv = clamp_llr((2.0 / sigma2) * y[v]);
        D->L[v] = Lv;
    }
    for (int m = 0; m < N; ++m) {
        for (int t = 0; t < DEG; ++t) {
            int v = D->check_to_vars[E(m, t)];
            D->q[E(m, t)] = D->L[v];
            D->r[E(m, t)] = 0.0;
        }
    }
}

static void make_aposteriori_llr(Decoder *D, double *L_post) {
    for (int v = 0; v < N; ++v) {
        double sum = D->L[v];
        for (int s = 0; s < DEG; ++s) {
            int m  = D->var_to_checks[v * DEG + s];
            int tp = D->var_pos_in_check[v * DEG + s];
            sum += D->r[E(m, tp)];
        }
        L_post[v] = clamp_llr(sum);
    }
}

static void llr_to_hard(const double *L_post, uint8_t *bits) {
    for (int v = 0; v < N; ++v) {
        bits[v] = (L_post[v] >= 0.0) ? 0u : 1u;
    }
}

static int syndrome_all_zero(Decoder *D, const uint8_t *hb) {
    for (int m = 0; m < N; ++m) {
        int acc = 0;
        for (int t = 0; t < DEG; ++t) {
            int v = D->check_to_vars[E(m, t)];
            acc ^= (hb[v] & 1);
        }
        if (acc != 0) return 0;
    }
    return 1;
}

// ============================= SPA / MSA =====================================
static void check_update_SPA(Decoder *D) {
    for (int m = 0; m < N; ++m) {
        double tvals[DEG], pref[DEG], suff[DEG];
        for (int t = 0; t < DEG; ++t) {
            double v = clamp_llr(D->q[E(m, t)]);
            tvals[t] = tanh(0.5 * v);
        }
        pref[0] = 1.0;
        for (int t = 1; t < DEG; ++t) {
            double p = pref[t - 1] * tvals[t - 1];
            pref[t] = clamp_val(p, -0.999999, 0.999999);
        }
        suff[DEG - 1] = 1.0;
        for (int t = DEG - 2; t >= 0; --t) {
            double p = suff[t + 1] * tvals[t + 1];
            suff[t] = clamp_val(p, -0.999999, 0.999999);
        }
        for (int t = 0; t < DEG; ++t) {
            double prod = clamp_val(pref[t] * suff[t], -0.999999, 0.999999);
            double rm = 2.0 * atanh(prod);
            D->r[E(m, t)] = clamp_llr(rm);
        }
    }
}

static void check_update_MSA(Decoder *D) {
    for (int m = 0; m < N; ++m) {
        double a[DEG];
        int sgn[DEG];
        for (int t = 0; t < DEG; ++t) {
            double v = clamp_llr(D->q[E(m, t)]);
            sgn[t] = (v >= 0.0) ? +1 : -1;
            a[t] = fabs(v);
        }
        for (int t = 0; t < DEG; ++t) {
            int sp = +1;
            double mn = 1e300;
            for (int u = 0; u < DEG; ++u) {
                if (u != t) {
                    sp *= sgn[u];
                    if (a[u] < mn) mn = a[u];
                }
            }
            double rm = (double)sp * mn;
            D->r[E(m, t)] = clamp_llr(rm);
        }
    }
}

static void variable_update(Decoder *D) {
    for (int v = 0; v < N; ++v) {
        double sum = D->L[v];
        double rb[DEG];
        for (int s = 0; s < DEG; ++s) {
            int m  = D->var_to_checks[v * DEG + s];
            int tp = D->var_pos_in_check[v * DEG + s];
            double rv = D->r[E(m, tp)];
            rb[s] = rv;
            sum += rv;
        }
        for (int s = 0; s < DEG; ++s) {
            int m  = D->var_to_checks[v * DEG + s];
            int tp = D->var_pos_in_check[v * DEG + s];
            D->q[E(m, tp)] = clamp_llr(sum - rb[s]);
        }
    }
}

static int decode_block_llr(Decoder *D, int max_iters, int algo, uint8_t *hard_bits, int *iters_used) {
    double *Lpost = (double*)malloc(sizeof(double) * N);
    if (!Lpost) { fprintf(stderr, "OOM Lpost\n"); exit(1); }
    for (int it = 1; it <= max_iters; ++it) {
        if (algo == 0) check_update_SPA(D);
        else           check_update_MSA(D);
        variable_update(D);
        make_aposteriori_llr(D, Lpost);
        llr_to_hard(Lpost, hard_bits);
        if (syndrome_all_zero(D, hard_bits)) {
            *iters_used = it;
            free(Lpost);
            return 1;
        }
    }
    *iters_used = max_iters;
    free(Lpost);
    return 0;
}

// ============================= Sim.txt =======================================
typedef struct {
    long long nBlocks;
    int       maxIters;
    double    ebn0_dB;
    long long seed;
    int       algo;
} SimConfig;

static int read_value_line(FILE *fp, char *line, size_t L) {
    if (!fgets(line, (int)L, fp)) return 0;
    for (char *p = line; *p; ++p) {
        if (*p == '%') { *p = '\0'; break; }
    }
    return 1;
}

static int load_Sim(SimConfig *cfg) {
    FILE *fp = fopen(FN_SIM, "r");
    if (!fp) {
        fprintf(stderr, "[ERR] open %s\n", FN_SIM);
        return 0;
    }
    char line[256];
    if (!read_value_line(fp, line, sizeof(line))) return 0; sscanf(line, "%lld", &cfg->nBlocks);
    if (!read_value_line(fp, line, sizeof(line))) return 0; sscanf(line, "%d",  &cfg->maxIters);
    if (!read_value_line(fp, line, sizeof(line))) return 0; sscanf(line, "%lf", &cfg->ebn0_dB);
    if (!read_value_line(fp, line, sizeof(line))) return 0; sscanf(line, "%lld", &cfg->seed);
    if (!read_value_line(fp, line, sizeof(line))) return 0; sscanf(line, "%d",  &cfg->algo);
    fclose(fp);
    return 1;
}

// ============================= main ==========================================
int main(void) {
    clock_t t0 = clock();

    // 讀參數
    SimConfig sim;
    if (!load_Sim(&sim)) {
        fprintf(stderr, "[ERR] 讀取 Sim.txt 失敗\n");
        return 1;
    }
    // 只支援 SPA(0)/MSA(1)；其他值一律改成 SPA
    if (sim.algo != 0 && sim.algo != 1) {
        fprintf(stderr, "[WRN] algo 非 0/1；以 SPA 執行\n");
        sim.algo = 0;
    }

    // CFG 起始列印（stdout）
    printf("CFG: blocks(max)=%lld, iters=%d, Eb/N0=%.2f dB, seed=%lld, algo=%s\n",
           sim.nBlocks, sim.maxIters, sim.ebn0_dB, sim.seed, (sim.algo==0?"SPA":"MSA"));
    if (STOP_BLOCK_ERRORS > 0ULL) {
        printf("CFG: early-stop when block_errors >= %llu\n",
               (unsigned long long)STOP_BLOCK_ERRORS);
    }
    if (STOP_BIT_ERRORS > 0ULL) {
        printf("CFG: early-stop when bit_errors   >= %llu\n",
               (unsigned long long)STOP_BIT_ERRORS);
    }
    fflush(stdout);

    // 讀取 G / H
    uint8_t **G = NULL;
    int *info_cols = NULL;
    if (!load_G(&G, &info_cols)) return 1;

    int *CTV = NULL, *VTC = NULL, *VPIC = NULL;
    if (!load_H(&CTV, &VTC, &VPIC)) return 1;

    Decoder D;
    if (!decoder_init(&D, CTV, VTC, VPIC)) return 1;

    // 通道參數與 seed（seed 為負數）
    double sigma2 = ebn0_db_to_sigma2(sim.ebn0_dB);
    double sigma  = sqrt(sigma2);
    long *idum = (long*)malloc(sizeof(long));
    *idum = (long)sim.seed;   // 由 Sim.txt 輸入，照老師規定為負數

    // 工作緩衝
    long long total_info_bits = sim.nBlocks * (long long)K;
    uint8_t *u_all = (uint8_t*)malloc(sizeof(uint8_t) * total_info_bits);
    uint8_t *u     = (uint8_t*)malloc(sizeof(uint8_t) * K);
    uint8_t *x     = (uint8_t*)malloc(sizeof(uint8_t) * N);
    double  *y     = (double *)malloc(sizeof(double ) * N);
    uint8_t *hard  = (uint8_t*)malloc(sizeof(uint8_t) * N);
    if (!u_all || !u || !x || !y || !hard) {
        fprintf(stderr, "[ERR] OOM main buffers\n");
        return 1;
    }

    make_information_bits(u_all, sim.nBlocks);

    // 統計量
    unsigned long long total_bits   = 0ULL;
    unsigned long long bit_errors   = 0ULL;
    unsigned long long total_blocks = 0ULL;
    unsigned long long block_errors = 0ULL;

    // 進度顯示
    long long step = sim.nBlocks / 10000;
    if (step < 1) step = 1;

    // 主迴圈（最多跑到 Sim.txt 設定的 blocks；達到門檻則提前停止）
    int early_stop = 0;
    const char *early_reason = NULL;
    long long early_at_block = -1;

    for (long long b = 0; b < sim.nBlocks; ++b) {
        // 取本 block 的資訊位
        for (int i = 0; i < K; ++i) {
            u[i] = u_all[b * (long long)K + i] & 1u;
        }

        // 編碼
        encode_block(G, u, x);

        // BPSK + AWGN（normal() 兩個一組）
        for (int j = 0; j < N; ) {
            double n1, n2;
            normal(&n1, &n2, sigma, idum);
            double s1 = bit_to_bpsk(x[j]);
            y[j] = s1 + n1;
            ++j;
            if (j < N) {
                double s2 = bit_to_bpsk(x[j]);
                y[j] = s2 + n2;
                ++j;
            }
        }

        // 解碼（只支援 SPA=0 / MSA=1）
        int iters_used = 0;
        initialize_llr_and_edges(&D, y, sigma2);
        (void)decode_block_llr(&D, sim.maxIters, sim.algo, hard, &iters_used);

        // 計錯（只看 systematic K 位）
        int err_in_block = 0;
        for (int i = 0; i < K; ++i) {
            uint8_t uhat = hard[ info_cols[i] ];
            uint8_t diff = (uint8_t)(uhat ^ u[i]);
            bit_errors += diff;
            err_in_block += diff;
        }
        total_bits += K;
        total_blocks += 1ULL;
        if (err_in_block > 0) block_errors += 1ULL;

        // 進度列（每 2% 或第一個 block 印一次；同時顯示目前錯的 block 數 / bit 數）
        if (((b + 1) % step) == 0 || b == 0) {
            double pct = 100.0 * (double)(b + 1) / (double)sim.nBlocks;
            printf("[INFO] %6.2f%%  (block %lld of %lld)  last iters=%d  block_errs=%llu  bit_errs=%llu\n",
                   pct,
                   (long long)(b + 1),
                   (long long)sim.nBlocks,
                   iters_used,
                   (unsigned long long)block_errors,
                   (unsigned long long)bit_errors);
            fflush(stdout);
        }

        // ---- 跨 block 早停 ----
        if (STOP_BLOCK_ERRORS > 0ULL && block_errors >= STOP_BLOCK_ERRORS) {
            early_stop = 1;
            early_reason = "block_errors threshold 100 reached";
            early_at_block = b + 1;
            break;
        }
        if (STOP_BIT_ERRORS > 0ULL && bit_errors >= STOP_BIT_ERRORS) {
            early_stop = 1;
            early_reason = "bit_errors threshold reached";
            early_at_block = b + 1;
            break;
        }
    }

    // 結果
    double BER  = (total_bits   > 0ULL) ? (double)bit_errors   / (double)total_bits   : 0.0;
    double BLER = (total_blocks > 0ULL) ? (double)block_errors / (double)total_blocks : 0.0;
    clock_t t1 = clock();
    double elapsed_s = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;

    if (early_stop) {
        printf("[INFO] early stop at block %lld (%s)\n",
               early_at_block,
               (early_reason ? early_reason : "threshold reached"));
    }

    printf("================ LDPC (1023,781) 模擬結果 ================\n");
    printf("Algo          : %s\n", (sim.algo == 0 ? "SPA" : "MSA"));
    printf("Blocks (cfg)  : %lld\n", sim.nBlocks);
    printf("Blocks (tried): %llu\n", (unsigned long long)total_blocks);
    printf("Max Iters     : %d\n", sim.maxIters);
    printf("Eb/N0 (dB)    : %.6f\n", sim.ebn0_dB);
    printf("sigma^2       : %.9f\n", sigma2);
    printf("Decoded bits  : %llu\n", (unsigned long long)total_bits);
    printf("Bit errors    : %llu\n", (unsigned long long)bit_errors);
    printf("BER           : %.9e\n", BER);
    printf("Block errors  : %llu\n", (unsigned long long)block_errors);
    printf("BLER          : %.9e\n", BLER);
    printf("Elapsed (sec) : %.3f\n", elapsed_s);
    printf("===========================================================\n");

    // 寫入 CSV（Algo 僅 SPA 或 MSA）
    ensure_csv_header("results_BER.csv");
    FILE *fcsv = fopen("results_BER.csv", "a");
    if (fcsv) {
        fprintf(fcsv, "%s,%.1f,%llu,%llu,%.9e,%llu,%.9e,%.3f\n",
                (sim.algo==0?"SPA":"MSA"),
                sim.ebn0_dB,
                (unsigned long long)total_bits,
                (unsigned long long)bit_errors,
                BER,
                (unsigned long long)total_blocks,
                BLER,
                elapsed_s);
        fclose(fcsv);
    }

    // 清理
    if (G) {
        for (int i = 0; i < K; ++i) if (G[i]) free(G[i]);
        free(G);
    }
    if (info_cols) free(info_cols);
    if (D.L) free(D.L);
    if (D.q) free(D.q);
    if (D.r) free(D.r);
    if (CTV)  free(CTV);
    if (VTC)  free(VTC);
    if (VPIC) free(VPIC);
    if (u_all) free(u_all);
    if (u) free(u);
    if (x) free(x);
    if (y) free(y);
    if (hard) free(hard);
    if (idum) free(idum);
    return 0;
}

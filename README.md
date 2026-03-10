# ldpc-code-c
C implementation of LDPC encoding and decoding

# LDPC (1023,781) C 語言實作

**LDPC (1023,781)** 的編解碼模擬，在 **AWGN 通道**下觀察 **BER / BLER** 表現。

## 專案內容

這個專案主要包含：
- LDPC (1023,781) 矩陣資料讀取
- 編碼 / 解碼流程實作
- AWGN 通道模擬
- BER / BLER 統計
- 結果圖輸出

## 我在這個專案中練到的能力

這個專案讓我比較扎實地練到幾個和韌體 / 底層開發很有關係的能力：

- **C 語言實作**
- **pointer 操作與陣列存取**
- **動態記憶體配置與釋放**
- **bitwise / XOR / parity 相關邏輯**
- **邊界檢查與除錯**
- **把數學演算法轉成可執行程式**
- **檢查模擬結果是否合理，並回頭修正程式流程**

對我來說，這個專案最有價值的地方，不只是把理論公式寫成程式，而是實際處理了很多工程上會遇到的問題，例如：

- pointer 與陣列存取
- 動態記憶體配置與釋放
- bit-level 邏輯處理（XOR / parity）
- 稀疏矩陣資料的讀取與使用
- 多層迭代流程的除錯
- 邊界檢查與結果驗證
---

## 專案結構

├─ src/                 # C source code
├─ data/                # LDPC matrix / input data
├─ results/             # BER / BLER 結果與圖
├─ report/              # 專案報告
├─ README.md

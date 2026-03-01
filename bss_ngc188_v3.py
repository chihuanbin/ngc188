import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from oc_tools_padova_edr3 import *

# ==========================================
# 1. 基础设置与数据读取
# ==========================================
plt.rcParams['xtick.minor.visible'], plt.rcParams['xtick.top'] = True, True
plt.rcParams['ytick.minor.visible'], plt.rcParams['ytick.right'] = True, True
plt.rcParams['xtick.direction'], plt.rcParams['ytick.direction'] = 'in', 'in'
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 16

# 读取数据 (请确保路径正确)
ucc_members = pd.read_csv('/Users/chihuanbin/Documents/ngc188/NGC_188_Filtered_Members_without_ruwe.csv')
ucc_members['bp_rp'] = ucc_members['phot_bp_mean_mag'] - ucc_members['phot_rp_mean_mag']
ucc_members['g_mag'] = ucc_members['phot_g_mean_mag']

# ==========================================
# 2. 生成等时线 (用于 BSS 参照和 RGB 轨迹)
# ==========================================
z = [9.8837, 1.806, 0.146, -0.03] # age(log), dist, Av, FeH
filters = ['Gmag', 'G_BPmag', 'G_RPmag']
refMag = 'Gmag'

load_mod_grid() 
grid_iso = get_iso_from_grid(z[0], (10. ** z[3]) * 0.0152, filters, refMag, nointerp=False)
fit_iso = make_obs_iso(filters, grid_iso, z[1], gaia_ext=True, Av=z[2])

cmd_iso = pd.DataFrame({
    'col': fit_iso['G_BPmag'] - fit_iso['G_RPmag'],
    'mag': fit_iso['Gmag']
})

# 确定 MSTO (主序转折点) 的颜色，这是 BSS 的核心参考
limt_color = cmd_iso['col'].min()
print(f"MSTO Color Limit: {limt_color:.4f}")

# ==========================================
# 3. 筛选逻辑
# ==========================================

# ------------------------------------------
# A. BSS 筛选 
# ------------------------------------------
# 定义：颜色比 MSTO 蓝，且亮度控制在 MSTO 上方到一定范围内
mask_bss = (ucc_members['bp_rp'] < limt_color) & \
           (ucc_members['g_mag'] < 15.7) & \
           (ucc_members['g_mag'] > 10.0) 
bss = ucc_members[mask_bss]

# ------------------------------------------
# B. RGB 科学筛选 (基于等时线轨迹)
# ------------------------------------------
# 1. 定义 RGB 演化起始点 (Base of RGB)
#    NGC 188 的次巨星分支(SGB)大约在 G=14.5 处结束并向上转折进入红巨星阶段
base_rgb_mag = 14.5

# 2. 提取等时线中的 RGB 部分建立模型
#    取亮度 < 14.5 且 颜色 > MSTO 的部分
iso_rgb_part = cmd_iso[(cmd_iso['mag'] < base_rgb_mag) & (cmd_iso['col'] > limt_color)].sort_values('mag')

# 3. 创建插值函数：输入亮度，返回理论颜色
#    f(g_mag) -> theoretical_bp_rp
func_rgb_track = interp1d(iso_rgb_part['mag'], iso_rgb_part['col'], 
                          kind='linear', fill_value='extrapolate')

# 4. 执行筛选
#    Step 1: 亮度必须亮于 RGB 底部
candidates = ucc_members[ucc_members['g_mag'] < base_rgb_mag].copy()

#    Step 2: 计算观测颜色与理论颜色的距离 (Residual)
candidates['theo_col'] = func_rgb_track(candidates['g_mag'])
candidates['col_diff'] = np.abs(candidates['bp_rp'] - candidates['theo_col'])

#    Step 3: 设定容差宽度 (Tolerance Width)
#    +/- 0.15 mag 可以包含光度误差、微小的红化差异以及红团簇(Red Clump)的宽度
rgb_width = 0.15

mask_rgb = (candidates['col_diff'] < rgb_width) & \
           (candidates['bp_rp'] > 0.8) # 额外的安全限制，防止选到极亮的主序星

rg = candidates[mask_rgb]

# 保存用于 K-S Test 的文件
bss.to_csv('ngc188_bss_final_no_ruwe.csv', index=False)
rg.to_csv('ngc188_rgb_precise_ruwe.csv', index=False)

print(f"筛选结果: BSS = {len(bss)} 颗, RGB = {len(rg)} 颗")

# ==========================================
# 4. 绘图验证
# ==========================================
fig, ax = plt.subplots(figsize=(8, 9))

# 1. 背景点
ax.scatter(ucc_members['bp_rp'], ucc_members['g_mag'], c='lightgrey', s=15, alpha=0.5, label='Members')

# 2. 绘制等时线 (黑色实线)
ax.plot(cmd_iso['col'], cmd_iso['mag'], c='k', linewidth=1.5, zorder=2, label='Isochrone')

# 3. 可视化 RGB 筛选通道 (红色虚线)
#    为了展示筛选的科学性，画出筛选边界
y_span = np.linspace(8, base_rgb_mag, 100)
x_center = func_rgb_track(y_span)
ax.plot(x_center - rgb_width, y_span, c='red', linestyle='--', linewidth=1, alpha=0.6)
ax.plot(x_center + rgb_width, y_span, c='red', linestyle='--', linewidth=1, alpha=0.6, label='RGB Selection Track')
#    画出 RGB 底部界限
ax.hlines(base_rgb_mag, 0.8, 1.8, colors='red', linestyles=':', linewidth=2, label='Base of RGB')

# 4. MSTO 界限 (蓝色虚线 - BSS 的参考)
ax.vlines(limt_color, 17, 9, colors='blue', linestyles='--', alpha=0.6, label='MSTO Limit')

# 5. 绘制选中的星
#    BSS: 蓝色五角星
ax.scatter(bss['bp_rp'], bss['g_mag'], marker='*', color='blue', s=120, edgecolors='k', zorder=10, label=f'BSS ({len(bss)})')
#    RGB: 红色圆点
ax.scatter(rg['bp_rp'], rg['g_mag'], marker='o', color='red', s=40, edgecolors='maroon', zorder=10, label=f'RGB ({len(rg)})')
# . ZAMS (如果有文件)
try:
    M67_ZAMS = pd.read_csv('/Users/chihuanbin/Documents/ngc188/BSS/ZAMS_M67.csv')
    d_DM = 11.28 - 9.614
    d_Ebp = (0.146 / (3.1 * 0.785)) - 0.054
    ax.plot(M67_ZAMS.bp_rp + d_Ebp, M67_ZAMS.Gmag + d_DM, c='gray', alpha=0.8, label='ZAMS')
except:
    print("未找到 ZAMS 文件，跳过绘制。")
# 6. 修饰
ax.invert_yaxis()
ax.set_xlim(0.1, 2.0)
ax.set_ylim(19, 8)
ax.set_xlabel('$G_{BP} - G_{RP}$')
ax.set_ylabel('$G$ magnitude')
ax.set_title('NGC 188: BSS vs Precise RGB Selection')

ax.legend(loc='lower left', fontsize=12, frameon=True, framealpha=0.9)

plt.tight_layout()
plt.savefig('ngc188_bss_vs_precise_rgb_no_ruwe.pdf')
plt.show()

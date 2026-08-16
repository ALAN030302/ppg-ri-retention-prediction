import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import os
from datetime import datetime
import matplotlib

# ==================== 用户可调整参数（请在此修改） ====================
input_file = r"C:\Users\姚钱磊\Desktop\补充实验预测\验证\预测结果\新建 Microsoft Excel 工作表.xlsx"
output_dir = "analysis_charts_zw"
include_ri = False

# ==================== 全局样式设置 ====================
global_font_size = 16
boxplot_font_size = 20
size_scatter_all = (10, 5)
size_scatter_individual = (4, 4)
size_error_dist = (7, 6)
size_violin = (8, 6)
size_residual = (10, 6)
size_metrics_bar = (5, 4)
size_boxplot = (12, 6)
size_trend = (9 , 5)
size_correlation = (8, 6)
size_rank = (12, 5)

scatter_alpha = 0.7
scatter_size = 50
marker_shape = 'o'
line_width = 2
grid_alpha = 0.3
bar_width = 0.15
dpi = 1200
transparent_bg = True

colors = {
    'rt_smrt_pred': '#012f48',
    'rt_M1_pred': '#7a0101',
    'rt_M2_pred': '#035830',
    'rt_M3_pred': '#669aba',
    'rt_actual': '#4c4c4c',
    'rti_M3_pred': '#be1420'
}

model_names = {
    'rt_smrt_pred': 'Literature Model (60k)',
    'rt_M1_pred': 'Non-Literature Model',
    'rt_M2_pred': 'Literature Condition Model',
    'rt_M3_pred': 'Literature RI Model'
}
# ================================================================

os.makedirs(output_dir, exist_ok=True)

# ==================== 关键修改：中英文分别设置字体 ====================
# 英文数字使用 Times New Roman，中文使用宋体 SimSun
plt.rcParams['font.family'] = ['Times New Roman', 'SimSun']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = global_font_size
plt.rcParams['axes.labelsize'] = global_font_size + 2
plt.rcParams['axes.titlesize'] = global_font_size + 4
plt.rcParams['legend.fontsize'] = global_font_size
plt.rcParams['xtick.labelsize'] = global_font_size
plt.rcParams['ytick.labelsize'] = global_font_size
plt.rcParams['lines.linewidth'] = line_width

def save_fig(fig, path, **kwargs):
    try:
        fig.tight_layout()
    except:
        fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    fig.savefig(path, dpi=dpi, transparent=transparent_bg, bbox_inches='tight')
    plt.close(fig)

def load_data(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"文件不存在: {file_path}")
    if file_path.endswith('.csv'):
        df = pd.read_csv(file_path)
    else:
        df = pd.read_excel(file_path)
    required_cols = ['smiles', 'rt', 'rt_smrt_pred', 'rt_M1_pred',
                     'rt_M2_pred', 'rt_M3_pred', 'rti_M3_pred']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"缺少列: {missing}")
    return df

def calculate_metrics(y_true, y_pred):
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    mask = ~np.isnan(y_true) & ~np.isnan(y_pred)
    y_true, y_pred = y_true[mask], y_pred[mask]
    if len(y_true) == 0:
        return {'R2': np.nan, 'MAE': np.nan, 'RMSE': np.nan, 'Pearson_r': np.nan}
    ss_res = np.sum((y_true - y_pred)**2)
    ss_tot = np.sum((y_true - np.mean(y_true))**2)
    r2 = 1 - ss_res/ss_tot if ss_tot != 0 else np.nan
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    if len(y_true) > 1 and np.std(y_true) > 0 and np.std(y_pred) > 0:
        pearson_r, _ = stats.pearsonr(y_true, y_pred)
    else:
        pearson_r = np.nan
    return {'R2': r2, 'MAE': mae, 'RMSE': rmse, 'Pearson_r': pearson_r}

def create_scatter_all_models(df, metrics, selected_models, save_path):
    fig, ax = plt.subplots(figsize=size_scatter_all)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    for model in selected_models:
        color = colors[model]
        label = model_names[model]
        r2 = metrics[model]['R2']
        label_text = f'{label} (R²={r2:.3f})' if not np.isnan(r2) else f'{label} (R²=NaN)'
        x = df['rt']; y = df[model]
        mask = ~x.isna() & ~y.isna()
        ax.scatter(x[mask], y[mask], alpha=scatter_alpha, s=scatter_size,
                   marker=marker_shape, color=color, label=label_text)
        if len(x[mask]) > 1 and np.std(y[mask]) > 0:
            z = np.polyfit(x[mask], y[mask], 1)
            p = np.poly1d(z)
            x_range = np.linspace(x[mask].min(), x[mask].max(), 100)
            ax.plot(x_range, p(x_range), '--', color=color, linewidth=line_width*0.7, alpha=0.7)
    min_val = min(df['rt'].min(), min([df[m].min() for m in selected_models]))
    max_val = max(df['rt'].max(), max([df[m].max() for m in selected_models]))
    ax.plot([min_val, max_val], [min_val, max_val], 'k:', alpha=0.3, label='理想线')
    ax.set_xlabel('实际保留时间')
    ax.set_ylabel('预测保留时间')
    ax.set_title('所有模型散点图')
    ax.grid(True, alpha=grid_alpha)
    ax.legend(loc='upper left', bbox_to_anchor=(1,1), frameon=False)
    save_fig(fig, save_path)

def create_individual_scatter(df, metrics, model_key, save_path):
    fig, ax = plt.subplots(figsize=size_scatter_individual)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    color = colors[model_key]
    name = model_names[model_key]
    x = df['rt']; y = df[model_key]
    mask = ~x.isna() & ~y.isna()
    ax.scatter(x[mask], y[mask], alpha=scatter_alpha, s=scatter_size,
               marker=marker_shape, color=color)
    if len(x[mask]) > 1 and np.std(y[mask]) > 0:
        z = np.polyfit(x[mask], y[mask], 1)
        p = np.poly1d(z)
        x_range = np.linspace(x[mask].min(), x[mask].max(), 100)
        ax.plot(x_range, p(x_range), 'r--', linewidth=line_width*1.5,
                label=f'y={z[0]:.3f}x+{z[1]:.3f}')
    min_val = min(x[mask].min(), y[mask].min())
    max_val = max(x[mask].max(), y[mask].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k:', alpha=0.3, label='理想线')
    r2 = metrics[model_key]['R2']
    mae = metrics[model_key]['MAE']
    text = f'R² = {r2:.3f}\nMAE = {mae:.2f}' if not np.isnan(r2) else 'R² = NaN'
    ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=global_font_size,
            verticalalignment='top', bbox=None)
    ax.set_xlabel('实际保留时间')
    ax.set_ylabel('预测保留时间')
    ax.set_title(f'{name} 的预测')
    ax.grid(True, alpha=grid_alpha)
    if 'regression_line' in locals():
        ax.legend(loc='lower right', frameon=False)
    save_fig(fig, save_path)

def create_error_distribution(df, selected_models, save_path):
    fig, ax = plt.subplots(figsize=size_error_dist)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    for model in selected_models:
        color = colors[model]
        label = model_names[model]
        err = df[model] - df['rt']
        mask = ~err.isna()
        ax.hist(err[mask], bins=30, alpha=0.5, color=color, label=label, density=True, edgecolor='black')
    ax.axvline(x=0, color='red', linestyle='--', linewidth=line_width, label='零误差')
    ax.set_xlabel('预测误差 (预测值 - 实际值)')
    ax.set_ylabel('密度')
    ax.set_title('误差分布')
    ax.grid(True, alpha=grid_alpha)
    ax.legend(frameon=False)
    save_fig(fig, save_path)

def create_absolute_error_violin(df, selected_models, save_path):
    fig, ax = plt.subplots(figsize=size_violin)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    data = []
    labels = []
    for model in selected_models:
        err = np.abs(df[model] - df['rt'])
        mask = ~err.isna()
        data.append(err[mask])
        labels.append(model_names[model])
    vp = ax.violinplot(data, showmeans=True, showmedians=True)
    for i, pc in enumerate(vp['bodies']):
        pc.set_facecolor(colors[selected_models[i]])
        pc.set_alpha(0.7)
    ax.set_xticks(range(1, len(labels)+1))
    ax.set_xticklabels(labels, rotation=15)
    ax.set_ylabel('绝对误差')
    ax.set_title('绝对误差分布 (小提琴图)')
    ax.grid(True, alpha=grid_alpha, axis='y')
    save_fig(fig, save_path)

def create_residual_plot(df, selected_models, save_path):
    fig, ax = plt.subplots(figsize=size_residual)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    for model in selected_models:
        color = colors[model]
        label = model_names[model]
        residual = df[model] - df['rt']
        mask = ~residual.isna()
        ax.scatter(df['rt'][mask], residual[mask], alpha=scatter_alpha*0.7,
                   s=scatter_size*0.7, color=color, label=label, marker=marker_shape)
    ax.axhline(y=0, color='red', linestyle='--', linewidth=line_width)
    ax.set_xlabel('实际保留时间')
    ax.set_ylabel('残差')
    ax.set_title('残差图')
    ax.grid(True, alpha=grid_alpha)
    ax.legend(loc='upper left', bbox_to_anchor=(1,1), frameon=False)
    save_fig(fig, save_path)

def create_metrics_bar_chart(metrics, selected_models, metric_name, save_path):
    fig, ax = plt.subplots(figsize=size_metrics_bar)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    values = [metrics[m][metric_name] for m in selected_models]
    labels = [model_names[m] for m in selected_models]
    x = np.arange(len(labels))
    bars = ax.bar(x, values, width=bar_width, color=[colors[m] for m in selected_models], alpha=0.7)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x()+bar.get_width()/2, val+0.01, f'{val:.3f}', ha='center', va='bottom', fontsize=global_font_size-1)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15)
    ax.set_ylabel(metric_name)
    ax.set_title(f'{metric_name} 比较')
    ax.grid(True, alpha=grid_alpha, axis='y')
    if metric_name == 'R2':
        ax.axhline(y=1.0, color='green', linestyle=':', alpha=0.5, label='完美')
        ax.axhline(y=0.0, color='red', linestyle=':', alpha=0.5, label='均值基线')
    save_fig(fig, save_path)

def create_boxplot(df, selected_models, save_path):
    fig, ax = plt.subplots(figsize=size_boxplot)
    if transparent_bg:
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    data = []
    means = []
    for model in selected_models:
        err = np.abs(df[model] - df['rt'])
        mask = ~err.isna()
        data.append(err[mask])
        means.append(np.mean(err[mask]))

    letters = [chr(ord('A') + i) for i in range(len(selected_models))]

    import matplotlib
    if int(matplotlib.__version__.split('.')[0]) >= 3 and int(matplotlib.__version__.split('.')[1]) >= 9:
        bp = ax.boxplot(data, tick_labels=letters, patch_artist=True)
    else:
        bp = ax.boxplot(data, labels=letters, patch_artist=True)

    for i, box in enumerate(bp['boxes']):
        box.set_facecolor(colors[selected_models[i]])
        box.set_alpha(0.7)

    legend_labels = [f'{letters[i]}: {model_names[model]} (均值={means[i]:.3f})'
                     for i, model in enumerate(selected_models)]
    import matplotlib.patches as mpatches
    patches = [mpatches.Patch(color=colors[model], label=label, alpha=0.7)
               for model, label in zip(selected_models, legend_labels)]
    ax.legend(handles=patches, loc='upper left', bbox_to_anchor=(1, 1),
              fontsize=boxplot_font_size - 1, frameon=False)

    ax.set_ylabel('绝对误差', fontsize=boxplot_font_size)
    ax.set_title('绝对误差箱线图', fontsize=boxplot_font_size + 2)
    ax.tick_params(axis='both', labelsize=boxplot_font_size)
    ax.set_xticklabels(letters, rotation=0, ha='center', fontsize=boxplot_font_size)

    ax.grid(True, alpha=grid_alpha, axis='y')
    fig.subplots_adjust(right=0.7, bottom=0.15)
    save_fig(fig, save_path)

def create_trend_comparison(df, selected_models, save_path):
    fig, ax = plt.subplots(figsize=size_trend)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    df_sorted = df.sort_values('rt').dropna(subset=['rt'] + selected_models)
    x = range(len(df_sorted))
    ax.plot(x, df_sorted['rt'], color=colors['rt_actual'], linewidth=line_width*1.5, label='实际值', zorder=5)
    for model in selected_models:
        ax.plot(x, df_sorted[model], '--', color=colors[model], alpha=0.7, linewidth=line_width, label=model_names[model])
    ax.set_xlabel('化合物 (按实际保留时间排序)')
    ax.set_ylabel('保留时间')
    ax.set_title('保留时间趋势比较')
    ax.grid(True, alpha=grid_alpha)
    ax.legend(loc='upper left', bbox_to_anchor=(1,1), frameon=False)
    save_fig(fig, save_path)

def create_correlation_matrix(df, selected_models, include_ri, save_path):
    cols = ['rt'] + selected_models
    if include_ri and 'rti_M3_pred' in df.columns:
        cols.append('rti_M3_pred')
    corr = df[cols].dropna().corr()
    rename = {'rt': '实际值'}
    rename.update(model_names)
    if include_ri:
        rename['rti_M3_pred'] = '保留指数'
    corr = corr.rename(index=rename, columns=rename)
    fig, ax = plt.subplots(figsize=size_correlation)
    if transparent_bg:
        fig.patch.set_alpha(0); ax.patch.set_alpha(0)
    im = ax.imshow(corr.values, cmap='coolwarm', vmin=-1, vmax=1)
    for i in range(len(corr)):
        for j in range(len(corr)):
            ax.text(j, i, f'{corr.iloc[i,j]:.3f}', ha='center', va='center',
                    color='white' if abs(corr.iloc[i,j])>0.5 else 'black', fontsize=global_font_size-1)
    ax.set_xticks(range(len(corr)))
    ax.set_yticks(range(len(corr)))
    ax.set_xticklabels(corr.columns, rotation=45, ha='right')
    ax.set_yticklabels(corr.index)
    ax.set_title('相关性矩阵')
    fig.colorbar(im, ax=ax)
    save_fig(fig, save_path)

def create_rank_comparison(df, selected_models, save_path):
    df_valid = df[['rt'] + selected_models].dropna()
    if len(df_valid) == 0:
        return
    ranks_list = []
    for _, row in df_valid.iterrows():
        actual = row['rt']
        errors = {m: abs(row[m] - actual) for m in selected_models}
        ranked = sorted(errors.items(), key=lambda x: x[1])
        ranks = {m: i+1 for i, (m, _) in enumerate(ranked)}
        ranks_list.append(ranks)
    avg_ranks = {m: np.mean([r[m] for r in ranks_list]) for m in selected_models}
    sorted_models = sorted(avg_ranks.items(), key=lambda x: x[1])
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=size_rank)
    if transparent_bg:
        fig.patch.set_alpha(0)
        for ax in (ax1, ax2):
            ax.patch.set_alpha(0)
    x = range(len(sorted_models))
    values = [v for _,v in sorted_models]
    labels = [model_names[m] for m,_ in sorted_models]
    bars = ax1.bar(x, values, width=0.6, color=[colors[m] for m,_ in sorted_models], alpha=0.7)
    for bar, val in zip(bars, values):
        ax1.text(bar.get_x()+bar.get_width()/2, val+0.05, f'{val:.2f}', ha='center', va='bottom')
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=15)
    ax1.set_ylabel('平均排名 (1=最佳)')
    ax1.set_title('模型平均性能排名')
    ax1.grid(True, alpha=grid_alpha, axis='y')
    n_models = len(selected_models)
    rank_mat = np.zeros((n_models, n_models))
    for r in ranks_list:
        for i, m in enumerate(selected_models):
            rank_mat[i, r[m]-1] += 1
    rank_mat = rank_mat / rank_mat.sum(axis=1, keepdims=True)
    im = ax2.imshow(rank_mat, cmap='YlOrRd', aspect='auto')
    for i in range(n_models):
        for j in range(n_models):
            if rank_mat[i,j] > 0:
                ax2.text(j, i, f'{rank_mat[i,j]:.2f}', ha='center', va='center',
                         color='white' if rank_mat[i,j]>0.5 else 'black', fontsize=global_font_size-2)
    ax2.set_xticks(range(n_models))
    ax2.set_yticks(range(n_models))
    ax2.set_xticklabels([f'排名 {i+1}' for i in range(n_models)])
    ax2.set_yticklabels([model_names[m] for m in selected_models])
    ax2.set_xlabel('排名位置')
    ax2.set_ylabel('预测模型')
    ax2.set_title('排名分布热图')
    fig.colorbar(im, ax=ax2, label='频率')
    save_fig(fig, save_path)

def main():
    print(f"输入文件: {input_file}")
    print(f"字体大小: {global_font_size}")
    print(f"输出目录: {output_dir}")
    df = load_data(input_file)
    print(f"成功加载 {len(df)} 个化合物。")
    selected_models = [c for c in model_names.keys() if c in df.columns]
    if not selected_models:
        raise ValueError("未找到有效的模型列。")
    metrics = {}
    for model in selected_models:
        metrics[model] = calculate_metrics(df['rt'], df[model])
    print("\n模型性能指标:")
    for model in selected_models:
        m = metrics[model]
        print(f"{model_names[model]}: R²={m['R2']:.4f}, MAE={m['MAE']:.4f}, RMSE={m['RMSE']:.4f}")

    print("\n绝对误差均值（各模型）:")
    for model in selected_models:
        err = np.abs(df[model] - df['rt'])
        mean_err = np.mean(err.dropna())
        print(f"{model_names[model]}: {mean_err:.4f}")

    print("\n正在生成图表...")
    create_scatter_all_models(df, metrics, selected_models, os.path.join(output_dir, '1_scatter_all_models.png'))
    for model in selected_models:
        create_individual_scatter(df, metrics, model, os.path.join(output_dir, f'2_scatter_{model}.png'))
    create_error_distribution(df, selected_models, os.path.join(output_dir, '3_error_distribution.png'))
    create_absolute_error_violin(df, selected_models, os.path.join(output_dir, '4_absolute_error_violin.png'))
    create_residual_plot(df, selected_models, os.path.join(output_dir, '5_residual_plot.png'))
    for metric_name in ['R2', 'MAE', 'RMSE', 'Pearson_r']:
        if all(not np.isnan(metrics[m][metric_name]) for m in selected_models):
            create_metrics_bar_chart(metrics, selected_models, metric_name,
                                     os.path.join(output_dir, f'6_bar_{metric_name}.png'))
    create_boxplot(df, selected_models, os.path.join(output_dir, '7_boxplot.png'))
    create_trend_comparison(df, selected_models, os.path.join(output_dir, '8_trend_comparison.png'))
    create_correlation_matrix(df, selected_models, include_ri, os.path.join(output_dir, '9_correlation_matrix.png'))
    create_rank_comparison(df, selected_models, os.path.join(output_dir, '10_rank_comparison.png'))
    print(f"\n所有图表已保存至 '{output_dir}' 目录。")

if __name__ == '__main__':
    main()
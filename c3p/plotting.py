import numpy as np
from matplotlib import pyplot as plt


def plot_scatter(df, x, y, title="Scatterplot", results_dir=None, log_x=False):
    """
    Create a scatter plot with optional logarithmic x-axis.

    Args:
        df: DataFrame containing the data
        x: Name of x-axis column
        y: Name of y-axis column
        title: Plot title
        results_dir: Optional directory to save the plot
        log_x: If True, x-axis will be logarithmic
    """
    correlation = np.corrcoef(df[x], df[y])[0, 1]
    x_values = np.log10(df[x]) if log_x else df[x]

    # Calculate correlation using transformed values if log scale
    #correlation = np.corrcoef(x_values, df[y])[0, 1]
    title = f"{title} (r = {correlation:.2f})"

    # Create scatter plot
    if log_x:
        plt.semilogx(df[x], df[y], 'o', markersize=5, alpha=0.5)
    else:
        plt.scatter(df[x], df[y], s=5, alpha=0.5)

    # Draw horizontal line at mean y
    plt.axhline(df[y].mean(), color='red', linestyle='dashed', linewidth=1)

    # Fit trend line using transformed x values
    z = np.polyfit(x_values, df[y], 1)
    p = np.poly1d(z)

    # Generate x points for trend line
    if log_x:
        x_trend = np.logspace(np.log10(df[x].min()), np.log10(df[x].max()), 100)
        plt.plot(x_trend, p(np.log10(x_trend)), "r--", alpha=0.8)
    else:
        plt.plot(df[x], p(df[x]), "r--", alpha=0.8)

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)
    plt.grid(True, alpha=0.3)

    if results_dir:
        scale_suffix = "-log" if log_x else ""
        plt.savefig(results_dir / f"scatterplot-{x}-{y}{scale_suffix}.png",
                    dpi=300,  # High resolution
                    bbox_inches='tight',  # Removes extra white spaces
                    pad_inches=0.1,  # Small padding around the figure
                    format='png',  # Format type
                    transparent=False,  # White background
                    facecolor='white',  # Figure face color
                    edgecolor='none',  # No edge color
                    )
    return plt
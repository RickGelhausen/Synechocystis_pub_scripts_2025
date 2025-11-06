"""
Contains scripts related to plotting metagene profiles.
Author: Rick Gelhausen
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots

def prepare_colors_and_lines(color_list):
    """
    Prepare a list of colors and a list of different line types for the plot.
    Ensures enough combinations for many read lengths.
    """
    colors = []
    lines = []
    line_types = ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]

    for line in line_types:
        colors.extend(color_list)
        lines.extend([line] * len(color_list))

    return colors, lines


def plot_metagene_from_dataframes(df_start, df_stop, chromosome, read_length_list, color_list=None):
    """
    Create metagene plot from dataframes with colorblind-friendly styling.
    """
    # Colorblind-friendly palette (works for deuteranopia, protanopia, and tritanopia)
    # Based on Paul Tol's colorblind-safe color schemes
    if color_list is None or color_list == []:
        color_list = [
            "#0173B2",  # Blue
            "#DE8F05",  # Orange
            "#029E73",  # Green
            "#CC78BC",  # Purple/Pink
            "#CA9161",  # Brown
            "#FBAFE4",  # Light pink
            "#949494",  # Gray
            "#ECE133",  # Yellow
            "#56B4E9",  # Sky blue
            "#D55E00"   # Vermillion
        ]

    # Prepare colors and line styles
    colors, lines = prepare_colors_and_lines(color_list)

    # Create subplots with better spacing
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Start Codon", "Stop Codon"),
        horizontal_spacing=0.12
    )

    # Calculate max_y for consistent y-axis
    max_y_start = 0
    max_y_stop = 0

    for read_length in read_length_list:
        df_rl_start = df_start[df_start['read_length'] == read_length]
        df_rl_stop = df_stop[df_stop['read_length'] == read_length]

        if not df_rl_start.empty:
            df_plot_start = df_rl_start.groupby('position')['count'].mean().reset_index()
            max_y_start = max(max_y_start, df_plot_start['count'].max())

        if not df_rl_stop.empty:
            df_plot_stop = df_rl_stop.groupby('position')['count'].mean().reset_index()
            max_y_stop = max(max_y_stop, df_plot_stop['count'].max())

    max_y = max(max_y_start, max_y_stop)
    max_y += max_y * 0.15  # Add 15% padding

    # Plot start codon
    color_idx = 0
    for read_length in read_length_list:
        df_rl = df_start[df_start['read_length'] == read_length]
        if df_rl.empty:
            continue

        df_plot = df_rl.groupby('position')['count'].mean().reset_index()

        fig.add_trace(
            go.Scatter(
                x=df_plot['position'],
                y=df_plot['count'],
                mode='lines+markers',  # Add markers for better distinction
                name=f"{read_length} nt",
                line=dict(
                    color=colors[color_idx],
                    dash=lines[color_idx],
                    width=2.5
                ),
                marker=dict(
                    size=4,
                    color=colors[color_idx],
                    symbol='circle'  # Different symbols can be added if needed
                ),
                legendgroup=str(read_length),
                showlegend=True,
                hovertemplate='<b>%{fullData.name}</b><br>' +
                             'Position: %{x}<br>' +
                             'Count: %{y:.1f}<br>' +
                             '<extra></extra>'
            ),
            row=1, col=1
        )
        color_idx += 1

    # Add vertical line at position 0 for start codon
    fig.add_shape(
        go.layout.Shape(
            type="line",
            x0=0, y0=0, x1=0, y1=max_y,
            line={"dash": "dash", "color": "rgba(0, 0, 0, 0.4)", "width": 2}
        ),
        row=1, col=1
    )

    # Plot stop codon
    color_idx = 0
    for read_length in read_length_list:
        df_rl = df_stop[df_stop['read_length'] == read_length]
        if df_rl.empty:
            continue

        df_plot = df_rl.groupby('position')['count'].mean().reset_index()

        fig.add_trace(
            go.Scatter(
                x=df_plot['position'],
                y=df_plot['count'],
                mode='lines+markers',
                name=f"{read_length} nt",
                line=dict(
                    color=colors[color_idx],
                    dash=lines[color_idx],
                    width=2.5
                ),
                marker=dict(
                    size=4,
                    color=colors[color_idx],
                    symbol='circle'
                ),
                legendgroup=str(read_length),
                showlegend=False,
                hovertemplate='<b>%{fullData.name}</b><br>' +
                             'Position: %{x}<br>' +
                             'Count: %{y:.1f}<br>' +
                             '<extra></extra>'
            ),
            row=1, col=2
        )
        color_idx += 1

    # Add vertical line at position 0 for stop codon
    fig.add_shape(
        go.layout.Shape(
            type="line",
            x0=0, y0=0, x1=0, y1=max_y,
            line={"dash": "dash", "color": "rgba(0, 0, 0, 0.4)", "width": 2}
        ),
        row=1, col=2
    )

    # Update axes with improved styling
    fig.update_xaxes(
        title_text="Distance to start codon (nt)",
        title_font={"size": 18, "family": "Arial, sans-serif", "color": "#2C3E50"},
        tickfont={"size": 14, "color": "#34495E"},
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(200, 200, 200, 0.3)',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='rgba(0, 0, 0, 0.2)',
        showline=True,
        linewidth=2,
        linecolor='rgba(0, 0, 0, 0.2)',
        mirror=True,
        row=1, col=1
    )

    fig.update_xaxes(
        title_text="Distance to stop codon (nt)",
        title_font={"size": 18, "family": "Arial, sans-serif", "color": "#2C3E50"},
        tickfont={"size": 14, "color": "#34495E"},
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(200, 200, 200, 0.3)',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='rgba(0, 0, 0, 0.2)',
        showline=True,
        linewidth=2,
        linecolor='rgba(0, 0, 0, 0.2)',
        mirror=True,
        row=1, col=2
    )

    fig.update_yaxes(
        title_text="Read Coverage",
        title_font={"size": 18, "family": "Arial, sans-serif", "color": "#2C3E50"},
        tickfont={"size": 14, "color": "#34495E"},
        range=[0, max_y],
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(200, 200, 200, 0.3)',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='rgba(0, 0, 0, 0.2)',
        showline=True,
        linewidth=2,
        linecolor='rgba(0, 0, 0, 0.2)',
        mirror=True,
        row=1, col=1
    )

    fig.update_yaxes(
        range=[0, max_y],
        tickfont={"size": 14, "color": "#34495E"},
        showgrid=True,
        gridwidth=1,
        gridcolor='rgba(200, 200, 200, 0.3)',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='rgba(0, 0, 0, 0.2)',
        showline=True,
        linewidth=2,
        linecolor='rgba(0, 0, 0, 0.2)',
        mirror=True,
        row=1, col=2
    )

    # Update subplot titles
    fig.update_annotations(
        font=dict(size=20, family="Arial, sans-serif", color="#2C3E50", weight='bold')
    )

    # Update overall layout
    fig.update_layout(
        title={
            "text": f"<b>{chromosome}</b>",
            "x": 0.5,
            "xanchor": "center",
            "font": {"size": 26, "family": "Arial, sans-serif", "color": "#1A252F"}
        },
        height=650,
        width=1400,
        showlegend=True,
        legend=dict(
            title=dict(text="Read Length", font=dict(size=16, family="Arial, sans-serif")),
            font=dict(size=14, family="Arial, sans-serif"),
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1,
            x=1.02,
            y=0.5,
            xanchor='left',
            yanchor='middle'
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=80, r=150, t=100, b=80),
        hovermode='x unified',
        font=dict(family="Arial, sans-serif")
    )

    return fig, max_y
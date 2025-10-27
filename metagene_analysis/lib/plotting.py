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
    Create metagene plot from dataframes.
    """
    # Set default colors if not provided
    if color_list is None or color_list == []:
        color_list = ["#4363d8", "#a9a9a9", "#ffe119", "#f58231", "#dcbeff",
                      "#000075", "#800000", "#e6194B", "#aaffc3", "#42d4f4"]

    # Prepare colors and line styles (creates many combinations)
    colors, lines = prepare_colors_and_lines(color_list)

    # Create subplots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Start Codon", "Stop Codon"),
        horizontal_spacing=0.1
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
    max_y += max_y * 0.1  # Add 10% padding

    # Plot start codon
    color_idx = 0
    for read_length in read_length_list:
        df_rl = df_start[df_start['read_length'] == read_length]
        if df_rl.empty:
            continue

        # Average across samples if multiple
        df_plot = df_rl.groupby('position')['count'].mean().reset_index()

        fig.add_trace(
            go.Scatter(
                x=df_plot['position'],
                y=df_plot['count'],
                mode='lines',
                name=f"{read_length} nt",
                line=dict(color=colors[color_idx], dash=lines[color_idx]),
                legendgroup=str(read_length),
                showlegend=True
            ),
            row=1, col=1
        )
        color_idx += 1

    # Add vertical line at position 0 for start codon
    fig.add_shape(
        go.layout.Shape(
            type="line",
            x0=0, y0=0, x1=0, y1=max_y,
            line={"dash": "dash", "color": "grey", "width": 2}
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
                mode='lines',
                name=f"{read_length} nt",
                line=dict(color=colors[color_idx], dash=lines[color_idx]),
                legendgroup=str(read_length),
                showlegend=False
            ),
            row=1, col=2
        )
        color_idx += 1

    # Add vertical line at position 0 for stop codon
    fig.add_shape(
        go.layout.Shape(
            type="line",
            x0=0, y0=0, x1=0, y1=max_y,
            line={"dash": "dash", "color": "grey", "width": 2}
        ),
        row=1, col=2
    )

    # Update axes
    fig.update_xaxes(
        title_text="Distance to start codon (nt)",
        title_font={"size": 20},
        row=1, col=1
    )
    fig.update_xaxes(
        title_text="Distance to stop codon (nt)",
        title_font={"size": 20},
        row=1, col=2
    )
    fig.update_yaxes(
        title_text="Read Coverage",
        title_font={"size": 20},
        range=[0, max_y],
        row=1, col=1
    )
    fig.update_yaxes(
        range=[0, max_y],
        row=1, col=2
    )

    # Update layout
    fig.update_annotations(font_size=24)
    fig.update_layout(
        title={
            "text": f"<b>{chromosome}</b>",
            "x": 0.5,
            "xanchor": "center",
            "font": {"size": 28}
        },
        height=600,
        showlegend=True
    )

    return fig, max_y
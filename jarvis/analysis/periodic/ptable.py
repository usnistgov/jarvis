"""Module for showing periodic table chemical trends."""

from bokeh.models import (
    ColumnDataSource,
    LinearColorMapper,
    LogColorMapper,
    ColorBar,
    BasicTicker,
)
from bokeh.plotting import figure
from bokeh.io import output_file, save, show
from bokeh.sampledata.periodic_table import elements
from bokeh.transform import dodge
from matplotlib.colors import Normalize, LogNorm, to_hex
from matplotlib.cm import plasma, ScalarMappable

# , inferno, magma, viridis, ScalarMappable

# Reuires bokeh installation https://docs.bokeh.org/en/latest/


def plot_ptable_trend(
    data_elements=["Rb", "S", "Se"],
    data_list=[10, 20, 30],
    input_file=None,
    output_html="ptable.html",
    bokeh_palette="Plasma256",
    cmap=plasma,
    log_scale=0,
    width=1050,
    alpha=0.65,
    cbar_height=520,
    cbar_font="14pt",
    save_plot=True,
):
    """
    Generate periodic table chemical trends.

    Either provide a file or list of data_elements, &data_list.
    Note that Bokeh already provided a periodic table.
    This module will take your data to color code them.
    See an example: https://www.nature.com/articles/s41598-019-45028-y
    Fig. 3
    Forked from https://github.com/arosen93/ptable_trends
    """
    output_file(output_html)
    # Define number of and groups
    period_label = ["1", "2", "3", "4", "5", "6", "7"]
    group_range = [str(x) for x in range(1, 19)]
    if input_file is not None:
        data_elements = []
        data_list = []

        f = open(input_file, "r")
        lines = f.read().splitlines()
        f.close()
        for i in lines:
            data_elements.append(i.split()[0])
            data_list.append(i.split()[1])

    data = [float(i) for i in data_list]

    if len(data) != len(data_elements):
        raise ValueError("Unequal number of atomic elements and data points")

    # lanthanides = [x.lower() for x in elements["symbol"][56:70].tolist()]
    # actinides = [x.lower() for x in elements["symbol"][88:102].tolist()]
    period_label.append("blank")
    period_label.append("La")
    period_label.append("Ac")

    count = 0
    for i in range(56, 70):
        elements.period[i] = "La"
        elements.group[i] = str(count + 4)
        count += 1

    count = 0
    for i in range(88, 102):
        elements.period[i] = "Ac"
        elements.group[i] = str(count + 4)
        count += 1

    # Define matplotlib and bokeh color map
    if log_scale == 0:
        color_mapper = LinearColorMapper(
            palette=bokeh_palette, low=min(data), high=max(data)
        )
        norm = Normalize(vmin=min(data), vmax=max(data))
    elif log_scale == 1:
        for i in range(len(data)):
            if data[i] < 0:
                raise ValueError(
                    "Entry for element "
                    + data_elements[i]
                    + " is negative but"
                    " log-scale is selected"
                )
        color_mapper = LogColorMapper(
            palette=bokeh_palette, low=min(data), high=max(data)
        )
        norm = LogNorm(vmin=min(data), vmax=max(data))
    color_scale = ScalarMappable(norm=norm, cmap=cmap).to_rgba(
        data, alpha=None
    )

    # Define color for blank entries
    blank_color = "#c4c4c4"
    color_list = []
    for i in range(len(elements)):
        color_list.append(blank_color)

    # Compare elements in dataset with elements in periodic table
    for i in range(len(data)):
        element_entry = elements.symbol[
            elements.symbol.str.lower() == data_elements[i].lower()
        ]
        if not element_entry.empty:
            element_index = element_entry.index[0]
        else:
            print("WARNING: Invalid chemical symbol: " + data_elements[i])
        if color_list[element_index] != blank_color:
            print("WARNING: Multiple entries for element " + data_elements[i])
        color_list[element_index] = to_hex(color_scale[i])

    # Define figure properties for visualizing data
    source = ColumnDataSource(
        data=dict(
            group=[str(x) for x in elements["group"]],
            period=[str(y) for y in elements["period"]],
            sym=elements["symbol"],
            atomic_number=elements["atomic number"],
            type_color=color_list,
        )
    )
    # Plot the periodic table
    p = figure(
        x_range=group_range, y_range=list(reversed(period_label)), tools="save"
    )
    p.plot_width = width
    p.outline_line_color = None
    p.toolbar_location = "above"
    p.rect(
        "group",
        "period",
        0.9,
        0.9,
        source=source,
        alpha=alpha,
        color="type_color",
    )
    p.axis.visible = False
    text_props = {
        "source": source,
        "angle": 0,
        "color": "black",
        "text_align": "left",
        "text_baseline": "middle",
    }
    x = dodge("group", -0.4, range=p.x_range)
    y = dodge("period", 0.3, range=p.y_range)
    p.text(
        x=x,
        y="period",
        text="sym",
        text_font_style="bold",
        text_font_size="15pt",
        **text_props
    )
    p.text(x=x, y=y, text="atomic_number", text_font_size="9pt", **text_props)
    color_bar = ColorBar(
        color_mapper=color_mapper,
        ticker=BasicTicker(desired_num_ticks=10),
        border_line_color=None,
        label_standoff=6,
        major_label_text_font_size=cbar_font,
        location=(0, 0),
        orientation="vertical",
        scale_alpha=alpha,
        width=8,
    )

    if cbar_height is not None:
        color_bar.height = cbar_height

    p.add_layout(color_bar, "right")
    p.grid.grid_line_color = None
    if save_plot:
        save(p)
    else:
        show(p)
    return p


"""
if __name__ == "__main__":
    x = plot_ptable_trend()
    from bokeh.io import export_svgs
    # For svg
    x.output_backend = "svg"
    export_svgs(x, filename="plot.svg")
"""

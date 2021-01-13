from dearpygui import core, simple
import uuid
import pytalises as pt
import numpy as np

core.set_main_window_title("pyTALISES")


def theme_callback(sender, data):
    core.set_theme(sender)


def initialize_wavefunction_config_values():
    core.add_value("num_ext_dim", 1)
    core.add_value("x_num_grid_points", 64)
    core.add_value("y_num_grid_points", 1)
    core.add_value("z_num_grid_points", 1)
    core.add_value("x_min", "-1")
    core.add_value("x_max", "+1")
    core.add_value("y_min", "-0")
    core.add_value("y_max", "+0")
    core.add_value("z_min", "-0")
    core.add_value("z_max", "+0")
    core.add_value("num_int_dim", 1)
    core.add_data("psi", ["0"])
    core.add_data("wavefunction_variables_names", [])
    core.add_data("wavefunction_variables_values", [])
    core.add_value("m", "1.054571817e-34")
    core.add_value("t0", "0.0")
    core.add_data("normalize_const", None)
    core.add_data("V", [])
    core.add_data("V_variables_names", [])
    core.add_data("V_variables_values", [])


def update_num_ext_dim(sender, data):
    core.set_value("num_ext_dim", core.get_value(sender) + 1)
    update_num_grid_points_input()
    update_spatial_ext_input()


def update_num_grid_points_input():
    num_ext_dim = core.get_value("num_ext_dim")
    core.delete_item(item="group_num_grid_points", children_only=True)
    if num_ext_dim >= 1:
        core.add_text("in x ", parent="group_num_grid_points")
        core.add_same_line(parent="group_num_grid_points")
        core.add_input_int(
            name="##x_num_grid_points",
            parent="group_num_grid_points",
            default_value=core.get_value("x_num_grid_points"),
            step=False,
            callback=update_num_grid_points,
        )
        if num_ext_dim == 1:
            core.set_value("y_num_grid_points", 1)
            core.set_value("z_num_grid_points", 1)
    if num_ext_dim >= 2:
        core.add_text("in y ", parent="group_num_grid_points")
        core.add_same_line(parent="group_num_grid_points")
        core.add_input_int(
            name="##y_num_grid_points",
            parent="group_num_grid_points",
            default_value=core.get_value("y_num_grid_points"),
            step=False,
            callback=update_num_grid_points,
        )
        if num_ext_dim == 2:
            core.set_value("z_num_grid_points", 1)
    if num_ext_dim >= 3:
        core.add_text("in z ", parent="group_num_grid_points")
        core.add_same_line(parent="group_num_grid_points")
        core.add_input_int(
            name="##z_num_grid_points",
            parent="group_num_grid_points",
            default_value=core.get_value("z_num_grid_points"),
            step=False,
            callback=update_num_grid_points,
        )


def update_num_grid_points(sender, data):
    if sender[2] == "x":
        core.set_value("x_num_grid_points", core.get_value(sender))
    if sender[2] == "y":
        core.set_value("y_num_grid_points", core.get_value(sender))
    if sender[2] == "z":
        core.set_value("z_num_grid_points", core.get_value(sender))


def update_spatial_ext_input():
    num_ext_dim = core.get_value("num_ext_dim")
    core.delete_item(item="group_spatial_ext", children_only=True)
    if num_ext_dim >= 1:
        core.add_text("min x", parent="group_spatial_ext")
        core.add_same_line(parent="group_spatial_ext")
        core.add_input_text(
            "##spatial_ext_min_x",
            parent="group_spatial_ext",
            callback=update_spatial_ext,
            default_value=str(core.get_value("x_min")),
            no_spaces=True,
            scientific=True,
        )
        core.add_text("max x", parent="group_spatial_ext")
        core.add_same_line(parent="group_spatial_ext")
        core.add_input_text(
            "##spatial_ext_max_x",
            parent="group_spatial_ext",
            callback=update_spatial_ext,
            default_value=str(core.get_value("x_max")),
            no_spaces=True,
            scientific=True,
        )
        if num_ext_dim == 1:
            core.set_value("y_min", "-0")
            core.set_value("y_max", "+0")
            core.set_value("z_min", "-0")
            core.set_value("z_max", "+0")
    if num_ext_dim >= 2:
        core.add_text("min y", parent="group_spatial_ext")
        core.add_same_line(parent="group_spatial_ext")
        core.add_input_text(
            "##spatial_ext_min_y",
            parent="group_spatial_ext",
            callback=update_spatial_ext,
            default_value=str(core.get_value("y_min")),
            no_spaces=True,
            scientific=True,
        )
        core.add_text("max y", parent="group_spatial_ext")
        core.add_same_line(parent="group_spatial_ext")
        core.add_input_text(
            "##spatial_ext_max_y",
            parent="group_spatial_ext",
            callback=update_spatial_ext,
            default_value=str(core.get_value("y_max")),
            no_spaces=True,
            scientific=True,
        )
        if num_ext_dim == 2:
            core.set_value("z_min", "-0")
            core.set_value("z_max", "+0")
    if num_ext_dim >= 3:
        core.add_text("min z", parent="group_spatial_ext")
        core.add_same_line(parent="group_spatial_ext")
        core.add_input_text(
            "##spatial_ext_min_z",
            parent="group_spatial_ext",
            callback=update_spatial_ext,
            default_value=str(core.get_value("z_min")),
            no_spaces=True,
            scientific=True,
        )
        core.add_text("max z", parent="group_spatial_ext")
        core.add_same_line(parent="group_spatial_ext")
        core.add_input_text(
            "##spatial_ext_max_z",
            parent="group_spatial_ext",
            callback=update_spatial_ext,
            default_value=str(core.get_value("z_max")),
            no_spaces=True,
            scientific=True,
        )


def update_spatial_ext(sender, data):
    if sender[-1] == "x":
        core.set_value("x_min", core.get_value("##spatial_ext_min_x"))
        core.set_value("x_max", core.get_value("##spatial_ext_max_x"))
    if sender[-1] == "y":
        core.set_value("y_min", core.get_value("##spatial_ext_min_y"))
        core.set_value("y_max", core.get_value("##spatial_ext_max_y"))
    if sender[-1] == "z":
        core.set_value("z_min", core.get_value("##spatial_ext_min_z"))
        core.set_value("z_max", core.get_value("##spatial_ext_max_z"))


def update_num_int_dim(sender, data):
    if core.get_value(sender) < 1:
        core.set_value(sender, 1)
        print("Number of internal dimensions must be greater than 0.")
        return
    core.set_value("num_int_dim", core.get_value(sender))
    update_wavefunction_amplitude_input()
    update_potential_matrix_input()


def update_wavefunction_amplitude_input():
    num_int_dim = core.get_value("num_int_dim")
    core.delete_item(item="group_wavefunction_amplitude", children_only=True)
    psi = core.get_data("psi")
    n_prev_elements = len(psi)
    n_new_elements = num_int_dim
    n_diff = n_new_elements - n_prev_elements
    if n_diff > 0:
        for i in range(n_diff):
            psi.append("0")

    if n_diff < 0:
        for i in range(-n_diff):
            psi.pop()

    for i in range(num_int_dim):
        core.add_input_text(
            f"##psi{i}",
            parent="group_wavefunction_amplitude",
            default_value=core.get_data("psi")[i],
            multiline=True,
            callback=update_wavefunction_amplitude,
            callback_data=i,
        )


def update_wavefunction_amplitude(sender, data):
    psi = core.get_data("psi")
    index = data
    psi[index] = core.get_value(sender)


def add_variable_input(sender, data):
    parent_group = data["parent_group"]
    variable_identifier = data["variable_identifier"]
    value_identifier = data["value_identifier"]
    variables_names = core.get_data(variable_identifier)
    variables_values = core.get_data(value_identifier)
    variable_id = len(variables_names)
    unique_name = uuid.uuid4()
    variables_names.append("")
    variables_values.append("")
    core.add_group(
        f"group_variable_{unique_name}",
        parent=parent_group,
        before=sender,
    )
    core.add_input_text(
        f"##variable_{unique_name}",
        width=60,
        default_value=core.get_data(variable_identifier)[variable_id],
        callback=update_variable_name,
        callback_data={"id": variable_id, "variable_identifier": variable_identifier},
    )
    core.add_same_line()
    core.add_text("=")
    core.add_same_line()
    core.add_input_text(
        f"##value_{unique_name}",
        width=60,
        default_value=core.get_data(value_identifier)[variable_id],
        no_spaces=True,
        scientific=True,
        callback=update_variable_value,
        callback_data={"id": variable_id, "value_identifier": value_identifier},
    )
    core.add_same_line()
    core.add_button(
        f"-##remove_variable_{unique_name}",
        callback=remove_variable_input,
        callback_data={
            "id": variable_id,
            "unique_name": unique_name,
            "variable_identifier": variable_identifier,
            "value_identifier": value_identifier,
        },
    )
    core.end(f"group_variable_{variable_id}")


def remove_variable_input(sender, data):
    variable_id = data["id"]
    unique_name = data["unique_name"]
    core.delete_item(f"group_variable_{unique_name}")
    variables_names = core.get_data(data["variable_identifier"])
    variables_values = core.get_data(data["value_identifier"])
    variables_names[variable_id] = ""
    variables_values[variable_id] = ""


def update_variable_name(sender, data):
    variable_id = data["id"]
    variables_names = core.get_data(data["variable_identifier"])
    variables_names[variable_id] = core.get_value(sender)


def update_variable_value(sender, data):
    variable_id = data["id"]
    variables_values = core.get_data(data["value_identifier"])
    variables_values[variable_id] = core.get_value(sender)


def generate_wavefunction_callback(sender, data):
    num_ext_dim = core.get_value("num_ext_dim")
    x_min = float(eval(core.get_value("x_min")))
    x_max = float(eval(core.get_value("x_max")))
    y_min = float(eval(core.get_value("y_min")))
    y_max = float(eval(core.get_value("y_max")))
    z_min = float(eval(core.get_value("z_min")))
    z_max = float(eval(core.get_value("z_max")))
    x_num_grid_points = core.get_value("x_num_grid_points")
    y_num_grid_points = core.get_value("y_num_grid_points")
    z_num_grid_points = core.get_value("z_num_grid_points")
    num_int_dim = core.get_value("num_int_dim")
    m = float(eval(core.get_value("m")))
    t0 = float(eval(core.get_value("t0")))
    normalize_const = core.get_value("normalize_const")
    print(normalize_const)
    if normalize_const is not None:
        normalize_const = float(eval(normalize_const))
    print(normalize_const)
    psi = core.get_data("psi")
    variables_names = [i for i in core.get_data("variables_names") if i != "" ]
    variables_values = [float(eval(i)) for i in core.get_data("variables_values") if i != ""]
    if variables_names is not None:
        variables = dict(zip(variables_names, variables_values))
    else:
        variables = {}
    wavefunction = pt.Wavefunction(
        psi=psi,
        number_of_grid_points=(x_num_grid_points, y_num_grid_points, z_num_grid_points),
        spatial_ext=[(x_min, x_max), (y_min, y_max), (z_min, z_max)],
        t0=t0,
        m=m,
        variables=variables,
        normalize_const=normalize_const,
    )
    core.add_data("wavefunction", wavefunction)
    core.delete_series(plot="##psi-plot", series="Amp")
    core.add_line_series(
        "##psi-plot", "Amp", list(wavefunction.r), list(np.abs(wavefunction.amp) ** 2)
    )


def update_potential_matrix_input():
    core.add_data(
        "V",
        ["0"]
        * (
            int(core.get_value("num_int_dim") / 2 * (core.get_value("num_int_dim") + 1))
        ),
    )
    core.delete_item("potential-matrix", children_only=True)
    core.add_group(str(uuid.uuid4()), parent="potential-matrix")
    for i in range(core.get_value("num_int_dim")):
        for j in range(i + 1):
            mat_element = f"V_{i}{j}"
            core.add_input_text(
                name=mat_element,
                label="",
                width=200,
                default_value=mat_element,
                callback=update_potential_matrix,
                callback_data=(i, j),
            )
            update_potential_matrix(mat_element, (i, j))
            if i != j:
                core.add_same_line()
    core.end()


def update_potential_matrix(sender, data):
    num_int_dim = core.get_value("num_int_dim")
    V = core.get_data("V")
    i, j = data
    index = i + j * num_int_dim - sum(range(j + 1))
    V[index] = core.get_value(sender)


with simple.window("pytalises-main-window"):
    pass

with simple.window("Settings", x_pos=0, y_pos=0, autosize=True):

    core.add_menu_bar("Menu")

    core.add_menu("Themes")
    core.add_menu_item("Dark", callback=theme_callback)
    core.add_menu_item("Light", callback=theme_callback)
    core.add_menu_item("Classic", callback=theme_callback)
    core.add_menu_item("Dark 2", callback=theme_callback)
    core.add_menu_item("Grey", callback=theme_callback)
    core.add_menu_item("Dark Grey", callback=theme_callback)
    core.add_menu_item("Cherry", callback=theme_callback)
    core.add_menu_item("Purple", callback=theme_callback)
    core.add_menu_item("Gold", callback=theme_callback)
    core.add_menu_item("Red", callback=theme_callback)
    core.end("Themes")

    core.add_menu("Tools")
    core.add_menu_item("Show Logger", callback=core.show_logger)
    core.add_menu_item("Show About", callback=simple.show_about)
    core.add_menu_item("Show Metrics", callback=simple.show_metrics)
    core.add_menu_item("Show Documentation", callback=simple.show_documentation)
    core.add_menu_item("Show Debug", callback=simple.show_debug)
    core.end("Tools")

    core.end("Menu")

    initialize_wavefunction_config_values()

    with simple.tab_bar(name="tab-menu"):
        with simple.tab(name="Hamiltonian"):

            # configuration of matrix elements
            core.add_text("Potential:")
            core.add_group("potential-matrix")
            update_potential_matrix_input()
            core.end("potential-matrix")

            # configuration for variables used in wavefunction
            core.add_text("Hamiltonian parameters:")
            core.add_group("group_hamiltonian_variables")
            core.add_button(
                "+##add_V_variable",
                callback=add_variable_input,
                callback_data={
                    "parent_group": "group_hamiltonian_variables",
                    "variable_identifier": "V_variables_names",
                    "value_identifier": "V_variables_values",
                },
            )
            core.end("group_hamiltonian_variables")

        with simple.tab(name="Wavefunction"):

            # configuration for the number of internal dimensions
            core.add_text("# of external dimensios:")
            core.add_radio_button(
                "buttons_num_ext_dim",
                items=[1, 2, 3],
                horizontal=True,
                callback=update_num_ext_dim,
            )

            # configuration for number of spatial grid points
            core.add_text("# of grid points:", tip="Should be a multiple of 2")
            core.add_group("group_num_grid_points")
            update_num_grid_points_input()
            core.end("group_num_grid_points")

            # configuration for spatial extent of wavefunction
            core.add_text("spatial extent:")
            core.add_group("group_spatial_ext")
            update_spatial_ext_input()
            core.end("group_spatial_ext")

            # configuration for number of internal dimensions
            core.add_text("# internal dimensions:")
            core.add_input_int(
                "##input_num_int_dim",
                default_value=core.get_value("num_int_dim"),
                callback=update_num_int_dim,
            )

            # configuration wave function ampltiudes
            core.add_text("wavefunction amplitudes:")
            core.add_group("group_wavefunction_amplitude")
            update_wavefunction_amplitude_input()
            core.end("group_wavefunction_amplitude")

            # configuration for variables used in wavefunction
            core.add_text("wavefunction parameters:")
            core.add_group("group_wavefunction_variables")
            core.add_button(
                "+##add_variable",
                callback=add_variable_input,
                callback_data={
                    "parent_group": "group_wavefunction_variables",
                    "variable_identifier": "wavefunction_variables_names",
                    "value_identifier": "wavefunction_variables_values",
                },
            )
            core.end("group_wavefunction_variables")

            # configuration for optional parameters
            # core.add_text("optional parameters:")
            core.add_collapsing_header("optional parameters")
            core.add_text("m = ")
            core.add_same_line()
            core.add_input_text(
                "kg",
                default_value=core.get_value("m"),
                no_spaces=True,
                scientific=True,
                callback=lambda sender, data: core.add_value(
                    "m", core.get_value(sender)
                ),
            )
            core.add_text("t0 =")
            core.add_same_line()
            core.add_input_text(
                "s",
                default_value=core.get_value("t0"),
                no_spaces=True,
                scientific=True,
                callback=lambda sender, data: core.add_value(
                    "t0", core.get_value(sender)
                ),
            )
            core.add_text("N0 =")
            core.add_same_line()
            core.add_input_text(
                "##normalization_constant",
                no_spaces=True,
                scientific=True,
                callback=lambda sender, data: core.add_value(
                    "normalize_const", core.get_value(sender)
                ),
            )
            core.end("optional parameters")
            core.add_spacing()

            # GENERATE WAVEFUNCTION
            core.add_button(
                "generate wavefunction",
                height=30,
                callback=generate_wavefunction_callback,
            )


with simple.window(
    name="psi-plot-windows",
    label="Wavefunction Amplitude",
    x_pos=300,
    y_pos=300,
    width=800,
    height=600,
):
    core.add_plot("##psi-plot")


core.start_dearpygui(primary_window="pytalises-main-window")

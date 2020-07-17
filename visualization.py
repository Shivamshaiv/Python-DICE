from economy_model import Economy
import plotly.graph_objects as go
import datetime


def make_datetime(arr, start_year=2005):
    x = []
    for i in range(len(arr)):
        x.append(datetime.datetime(year=start_year + (5 * i), month=1, day=1))
    return x


def plot_attr(
    arr, title="Evolution of the factor in the model", desc="Model Parameter"
):
    model = Economy()
    fig = go.Figure(
        data=go.Scatter(x=make_datetime(arr, start_year=model.start_year), y=arr),
        layout=go.Layout(title=title,),
    )

    fig.update_xaxes(title_text="Time in Years")
    fig.update_yaxes(title_text=desc)

    fig.update_layout(template="none")
    # fig = go.Figure(data=[go.Scatter(x=make_datetime(model.T_at), y=model.T_at)])
    # Use datetime objects to set xaxis range
    fig.update_layout(
        xaxis_range=[datetime.datetime(2001, 1, 1), make_datetime(arr)[-1]]
    )
    fig.show()


def plot_temp_change(decades=10, tipping_damage=False, carbon_tax=(0, 0, 0)):
    model = Economy()
    model.loop(t=decades, tipping_damage=tipping_damage, carbon_tax=carbon_tax)
    fig = go.Figure(
        data=go.Scatter(
            x=make_datetime(model.T_at, start_year=model.start_year), y=model.T_at
        ),
        layout=go.Layout(title="Average Temperature Rise",),
    )

    fig.update_xaxes(title_text="Time in Years")
    fig.update_yaxes(title_text="Rise in Average Temperature since 1900")

    fig.update_layout(template="none")
    # fig = go.Figure(data=[go.Scatter(x=make_datetime(model.T_at), y=model.T_at)])
    # Use datetime objects to set xaxis range
    fig.update_layout(
        xaxis_range=[datetime.datetime(2001, 1, 1), make_datetime(model.T_at)[-1]]
    )
    fig.show()


def plot_multiple_attr(main_arr):
    """
    Input an array in the format of [(arr_1,name_1),(arr_2,name_2),(arr_3,name_3),...,(arr_n,name_n)]
    Please ensure the sizes of all the arrays are same
    """
    fig = go.Figure()
    for tup in main_arr:
        arr, name = tup
        this_attr = go.Scatter(x=make_datetime(arr), y=arr, name=str(name))
        fig.add_trace(this_attr)
    fig.update_xaxes(title_text="Time in Years")
    fig.update_layout(template="none")
    fig.show()


def give_temperature(decades, tipping_damage=False, carbon_tax=(0, 0, 0)):
    model = Economy()
    model.loop(decades, tipping_damage=tipping_damage, carbon_tax=carbon_tax)
    return model.T_at


def plot_impact_of_carbon_tax(
    decades, tipping_damage=False, tax_arr=[], show_default=True
):
    """
    To see the impact of diffrent climate policies and represent them in one plot.
    Input : decades,tipping_damage,tax_arr
    Here tax_arr is an array of tuples with three values each : the amount of tax per topn of CO2 in 2050,2100
    and 2150.
    """
    fig = go.Figure()
    if show_default:
        default_arr = give_temperature(decades, tipping_damage=tipping_damage)
        fig.add_trace(
            go.Scatter(
                x=make_datetime(default_arr), y=default_arr, name="No Carbon Tax"
            )
        )

    # In the loop for the carbon tax senarios
    for carbon_tax_tup in tax_arr:
        arr = give_temperature(
            decades, tipping_damage=tipping_damage, carbon_tax=carbon_tax_tup
        )
        fig.add_trace(go.Scatter(x=make_datetime(arr), y=arr, name=str(carbon_tax_tup)))

    fig.update_xaxes(title_text="Time in Years")
    fig.update_layout(template="none")
    fig.show()

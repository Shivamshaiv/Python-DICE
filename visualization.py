from economy_model import Economy
import plotly.graph_objects as go
import datetime

def make_datetime(arr,start_year = 2005):
    x = []
    for i in range(len(arr)):
        x.append(datetime.datetime(year=start_year+(10*i),month=1,day=1))
    return x

def plot_attr(arr,title="Evolution of the factor in the model",desc="Model Parameter"):
    model = Economy()
    fig = go.Figure(
    data=go.Scatter(x=make_datetime(arr,start_year=model.start_year), y=arr),
    layout=go.Layout(
        title=title,
    ))

    fig.update_xaxes(title_text='Time in Years')
    fig.update_yaxes(title_text=desc)

    fig.update_layout(template="none")
    #fig = go.Figure(data=[go.Scatter(x=make_datetime(model.T_at), y=model.T_at)])
    # Use datetime objects to set xaxis range
    fig.update_layout(xaxis_range=[datetime.datetime(2001, 1, 1),make_datetime(arr)[-1]])
    fig.show()

def plot_temp_change(decades = 10,tipping_damage = False):
    model = Economy()
    model.loop(t = decades,tipping_damage = tipping_damage)
    fig = go.Figure(
    data=go.Scatter(x=make_datetime(model.T_at,start_year=model.start_year), y=model.T_at),
    layout=go.Layout(
        title="Average Temperature Rise",
    ))

    fig.update_xaxes(title_text='Time in Years')
    fig.update_yaxes(title_text='Rise in Average Temperature since 1900')

    fig.update_layout(template="none")
    #fig = go.Figure(data=[go.Scatter(x=make_datetime(model.T_at), y=model.T_at)])
    # Use datetime objects to set xaxis range
    fig.update_layout(xaxis_range=[datetime.datetime(2001, 1, 1),make_datetime(model.T_at)[-1]])
    fig.show()

def plot_multiple_attr(main_arr):
    '''
    Input an array in the format of [(arr_1,name_1),(arr_2,name_2),(arr_3,name_3),...,(arr_n,name_n)]
    Please ensure the sizes of all the arrays are same
    '''
    fig = go.Figure()
    for tup in main_arr:
        arr,name = tup
        this_attr = go.Scatter(x=make_datetime(arr), y=arr,name=str(name))
        fig.add_trace(this_attr)
    fig.update_xaxes(title_text='Time in Years')
    fig.update_layout(template="none")
    fig.show()

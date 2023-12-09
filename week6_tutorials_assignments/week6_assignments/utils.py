import requests
from typing import List, Union
from time import sleep
import simplejson

from evidently.ui.dashboards import ReportFilter, DashboardPanelTestSuiteCounter, CounterAgg, DashboardPanelPlot, PanelValue, PlotType, DashboardConfig
from evidently.renderers.html_widgets import WidgetSize
from evidently import metrics
from evidently.ui.remote import RemoteWorkspace, RemoteProject
from evidently.ui.workspace import Workspace, Project

def send_requests(model_name: str, input: List, count: int):
    """
    Send some requests to an inference service

    Args:
        model_name (str): Name of the inference service
        input: Features based on which the inference service make predictions
        count: Number of requests
    """
    request_data = {
        'parameters': {'content_type': 'pd'}, 
        'inputs': input
    }
    headers = {}
    kserve_gateway_url = "http://kserve-gateway.local:30200"
    
    headers["Host"] = f"{model_name}.kserve-inference.example.com"
    headers["content-type"]="application/json"

    url = f"{kserve_gateway_url}/v2/models/{model_name}/infer"
    if (count == 1):
        res = requests.post(url, data=simplejson.dumps(request_data, ignore_nan=True), headers=headers)
        return res
        
    for i in range(count):
        requests.post(url, data=simplejson.dumps(request_data, ignore_nan=True), headers=headers)
        sleep(0.5)
        print(f"{i+1} requests have been sent", flush=True)


def init_evidently_project(workspace: Union[Workspace, RemoteWorkspace], project_name: str) -> Union[Project,RemoteProject]:
    """
    Create a Project to a Workspace
    Args:
        workspace: An Evidently Workspace
        project_name: Name of the Project
    """
    # Delete any projects whose name is the given project_name to avoid duplicated projects
    for project in workspace.search_project(project_name=project_name):
        workspace.delete_project(project_id=project.id)

    # Create a project at Evidently
    project = workspace.create_project(name=project_name)

    # Create a dashboard
    project.dashboard = DashboardConfig(name=project_name, panels=[])

    project.dashboard.add_panel(
        DashboardPanelTestSuiteCounter(
            title="MAE",
            agg=CounterAgg.LAST
        ),
    )
    project.dashboard.add_panel(
         DashboardPanelPlot(
                title="MAE",
                filter=ReportFilter(metadata_values={}, tag_values=[]),
                values=[
                    PanelValue(
                        metric_id="RegressionQualityMetric",
                        field_path=metrics.RegressionQualityMetric.fields.current.mean_abs_error,
                        legend="MAE",
                    ),
                ],
                plot_type=PlotType.LINE,
                size=WidgetSize.FULL,
            )
    )

    project.save()
    return project

# if __name__ == "__main__":
#     arr = [{'name': 'yr_built', 'shape': [10, 1], 'datatype': 'FP64', 'data': [1989.0, 1976.0, 1968.0, 1997.0, 1950.0, 2014.0, 1972.0, 1984.0, 1989.0, 2006.0]}, {'name': 'bedrooms', 'shape': [10, 1], 'datatype': 'FP64', 'data': [4.0, 4.0, 5.0, 4.0, 4.0, 5.0, 3.0, 3.0, 3.0, 1.0]}, {'name': 'postcode', 'shape': [10, 1], 'datatype': 'FP64', 'data': [707747.0018370263, 597736.9217030165, 597736.9217030165, 474632.3790458373, 515187.232483403, 468060.0543995248, 573737.7786048551, 474632.3790458373, 554723.8083336328, 554723.8083336328]}, {'name': 'area', 'shape': [10, 1], 'datatype': 'FP64', 'data': [726350.8929052162, 573567.1620235543, 675932.4355670998, 481874.4106608018, 512186.3536409238, 467555.6662013869, 564951.1784201534, 462591.44565758924, 566238.286454928, 564169.3413942211]}, {'name': 'bathrooms', 'shape': [10, 1], 'datatype': 'FP64', 'data': [3.25, 2.0, 2.75, 2.0, 2.0, 3.0, 2.0, 2.0, 2.5, 1.25]}, {'name': 'condition', 'shape': [10, 1], 'datatype': 'FP64', 'data': [3.0, 4.0, 4.0, 3.0, 3.0, 2.0, 4.0, 3.0, 3.0, np.nan]}, {'name': 'grade', 'shape': [10, 1], 'datatype': 'FP64', 'data': [11.0, 9.0, 9.0, 9.0, 8.0, 9.0, 8.0, 8.0, 8.0, 7.0]}, {'name': 'sqft_living', 'shape': [10, 1], 'datatype': 'FP64', 'data': [4200.0, 3085.0, 3050.0, 3345.0, 2195.0, 3045.0, 1805.0, 2115.0, 1940.0, 930.0]}, {'name': 'sqft_lot', 'shape': [10, 1], 'datatype': 'FP64', 'data': [18729.0, 20616.0, 13079.0, 34672.0, 8057.0, 8972.0, 8216.0, 38640.0, 6751.0, 2235.0]}, {'name': 'sqft_basement', 'shape': [10, 1], 'datatype': 'FP64', 'data': [110.0, 1070.0, 100.0, 0.0, 810.0, 20.0, 150.0, 90.0, 10.0, 30.0]}, {'name': 'sqft_living15', 'shape': [10, 1], 'datatype': 'FP64', 'data': [4097.0, 3052.0, 3033.0, 3233.0, 2023.0, 2930.0, 2036.0, 2287.0, 2054.0, 1305.0]}, {'name': 'sqft_lot15', 'shape': [10, 1], 'datatype': 'FP64', 'data': [18395.0, 19540.0, 12360.0, 23993.0, 7454.0, 8993.0, 8061.714182174384, 36053.0, 6638.0, 2827.0]}, {'name': 'waterfront', 'shape': [10, 1], 'datatype': 'FP64', 'data': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {'name': 'view', 'shape': [10, 1], 'datatype': 'FP64', 'data': [0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}, {'name': 'distance', 'shape': [10, 1], 'datatype': 'FP64', 'data': [3.41, 11.11, 19.3, 18.91, 13.43, 11.87, 9.68, 18.05, 11.37, 8.89]}, {'name': 'year', 'shape': [10, 1], 'datatype': 'INT64', 'data': [2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018, 2018]}]
#     send_requests(model_name="bike-price", input=arr, count=1)





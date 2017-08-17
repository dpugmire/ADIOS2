/*
 * VisVTKm.cpp
 *
 *  Created on: Aug 16, 2017
 *      Author: wfg
 */

#include "VisVTKm.h"

#ifndef VTKM_DEVICE_ADAPTER
#define VTKM_DEVICE_ADAPTER VTKM_DEVICE_ADAPTER_SERIAL
#endif

#include <iostream>
#include <vtkm/Math.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/filter/MarchingCubes.h>

#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Color.h>
#include <vtkm/rendering/Mapper.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/TextAnnotationScreen.h>
#include <vtkm/rendering/View1D.h>
#include <vtkm/rendering/View2D.h>
#include <vtkm/rendering/View3D.h>
namespace adios2
{
/*
//template <>
void SetCamera<vtkm::rendering::View3D>(vtkm::rendering::Camera& camera,
                                        const vtkm::Bounds& coordBounds)
{
  camera = vtkm::rendering::Camera();
  camera.ResetToBounds(coordBounds);
  camera.Azimuth(static_cast<vtkm::Float32>(45.0));
  camera.Elevation(static_cast<vtkm::Float32>(45.0));
}


template <typename MapperType, typename CanvasType, typename ViewType>
void Render(ViewType& view, const std::string& outputFile)
{
  view.Initialize();
  view.Paint();
  view.SaveAs(outputFile);
}
*/

void Render(const vtkm::cont::DataSet& ds,
            const std::string& fieldNm,
            const vtkm::rendering::ColorTable& colorTable,
            const std::string& outputFile)
{
    vtkm::rendering::MapperRayTracer mapper;
    vtkm::rendering::CanvasRayTracer canvas(512, 512);
  canvas.SetBackgroundColor(vtkm::rendering::Color::white);
  vtkm::rendering::Scene scene;

  scene.AddActor(vtkm::rendering::Actor(
    ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), colorTable));
  vtkm::rendering::Camera camera;
  camera = vtkm::rendering::Camera();
  camera.ResetToBounds(ds.GetCoordinateSystem().GetBounds());
  camera.Azimuth(static_cast<vtkm::Float32>(45.0));
  camera.Elevation(static_cast<vtkm::Float32>(45.0));
  
  vtkm::rendering::View3D view(scene, mapper, canvas, camera, vtkm::rendering::Color(0.2f, 0.2f, 0.2f, 1.0f));

  //Render<MapperType, CanvasType, ViewType>(view, outputFile);
  view.Initialize();
  view.Paint();
  view.SaveAs(outputFile);
}


bool VisVTKm::RenderAllVariables()
{
    for (auto &visVariable : m_VisVariables)
    {
        auto &var = visVariable.VisVariable;
        const void *buff = visVariable.Data;
        std::cout<<"BUFF Size: "<<visVariable.Size<<std::endl;

        std::cout << "Variable name " << var.m_Name << std::endl;
        std::cout<<"SHAPE: " << var.m_Shape.size() << " : " << var.m_Shape[0] << std::endl;
        std::cout<<"START: " << var.m_Start.size() << " : " << var.m_Start[0] << std::endl;
        std::cout<<"COUNT: " << var.m_Count.size() << " : " << var.m_Count[0] << std::endl;
        
        // Create the dataset from the variables        
        vtkm::Id3 dims(var.m_Shape[0], var.m_Shape[1], var.m_Shape[2]);
        vtkm::cont::DataSetBuilderUniform dsb;
        vtkm::cont::DataSet ds = dsb.Create(dims);
        
        // Add field to ds
        // Get the actual variable data to create the field
        const float *varBuff = (const float *)buff;
        vtkm::Id numPoints = dims[0]*dims[1]*dims[2];
        
        vtkm::cont::DataSetFieldAdd dsf;
        dsf.AddPointField(ds, var.m_Name, varBuff, numPoints);
        //ds.PrintSummary(std::cout);
    
        for (auto &transform : var.m_TransformsInfo)
        {
            // transform parameters
            for (auto &param : transform.Operator.m_Parameters)
            {
                if(param.first == "iso")
                {
                    // How to handle multiple iso values? Need to change where executed
                    vtkm::filter::MarchingCubes filter;
                    filter.SetGenerateNormals(true);
                    filter.SetMergeDuplicatePoints(false);
                    filter.SetIsoValue(0, stof(param.second));
                    
                    vtkm::filter::ResultDataSet result = filter.Execute(ds, ds.GetField(var.m_Name));
                    vtkm::cont::DataSet& outputData = result.GetDataSet();
                    std::cout<<"***************************************************************"<<std::endl;
                    std::cout<<"***************************************************************"<<std::endl;
                    outputData.PrintSummary(std::cout);
                    std::cout<<"***************************************************************"<<std::endl;
                    std::cout<<"***************************************************************"<<std::endl;

                    vtkm::Id numIsoPts = outputData.GetCellSet(0).GetNumberOfPoints();
                    std::vector<float> isoVar(numIsoPts, stof(param.second));
                    dsf.AddPointField(outputData, var.m_Name, isoVar);

                    //Now, render the Mr. Isosurface
                    Render(outputData, var.m_Name, vtkm::rendering::ColorTable("thermal"), "mr_iso.pnm");
                    Render(ds, var.m_Name, vtkm::rendering::ColorTable("thermal"), "mr_data.pnm");
                    

                }
                std::cout<<param.first<<" "<<param.second<<std::endl;
            }
            
            std::cout << __LINE__ << std::endl;
            for (auto &parameter : transform.Parameters)
            {
                const std::string key(parameter.first);
                const std::string value(parameter.second);

                std::cout << parameter.first << "  " << parameter.second << std::endl;
                if (key == "X1")
                {
                    auto value = parameter.second;
                    std::cout << __LINE__ << std::endl;
                    std::cout << "Meow" << std::endl;
                    /// CAll VTKm magic
                }
            }
        }
    }
    
    return true;
}
}

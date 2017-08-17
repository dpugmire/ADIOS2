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

namespace adios2
{

bool VisVTKm::RenderAllVariables()
{
    for (auto &visVariable : m_VisVariables)
    {
        auto &var = visVariable.VisVariable;
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
        std::vector<float> field;
        for ( int i = 0; i < dims[0]*dims[1]*dims[2]; i++)
        {
            field.push_back(i+.25f);
        }
        
        vtkm::cont::DataSetFieldAdd dsf;
        dsf.AddPointField(ds, var.m_Name, field);

        ds.PrintSummary(std::cout);
    
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
                    
                    /*
                    vtkm::filter::Result result = filter.Execute(ds, ds.GetField(var.m_Name));
                    filter.MapFieldOntoOutput(result, ds.GetField(var.m_Name));

                    //need to extract vertices, normals, and scalars
                    vtkm::cont::DataSet& outputData = result.GetDataSet();
                    outputData.PrintSummary(std::cout);
                    */
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

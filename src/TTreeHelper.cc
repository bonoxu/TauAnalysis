#include <iostream>

#include <TTreeHelper.h>

//------------------------------------------------------------------------------------------------------------------------------------------

TTreeHelper::TTreeHelper(const std::string &tFileName, const std::string &tFileOption, const std::string &tTreeName, const std::string &tTreeID)
{
    m_pTFile  = new TFile(tFileName.c_str(), tFileOption.c_str());
    m_pTTree = new TTree(tTreeName.c_str(), tTreeID.c_str());
}

//------------------------------------------------------------------------------------------------------------------------------------------

TTreeHelper::~TTreeHelper()
{
    m_pTFile->Close();
    delete m_pTFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::FillNode()
{
    m_pTTree->Fill();
    this->Clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::WriteTTreeToTFile()
{
    m_pTFile->cd();
    m_pTTree->Write();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::Clear()
{
    for (StrIntMap::iterator iter = m_strIntMap.begin(), iterEnd = m_strIntMap.end(); iter != iterEnd; ++iter)
    {
        iter->second = 0;
    }
    for (StrDoubleMap::iterator iter = m_strDoubleMap.begin(), iterEnd = m_strDoubleMap.end(); iter != iterEnd; ++iter)
    {
        iter->second = 0.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::InitialiseIntVar(const std::string &varName)
{
    if (!m_strIntMap.insert(StrIntMap::value_type(varName, 0)).second)
        std::cout<<"TTreeHelper::InitialiseIntVar error m_strIntMap unable to insert "<<varName<<std::endl;
    m_pTTree->Branch(varName.c_str(), &m_strIntMap.at(varName), "varName/I");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::InitialiseIntVars(const StrSet &varNames)
{
    for (StrSet::const_iterator iter = varNames.begin(), iterEnd = varNames.end(); iter != iterEnd; ++iter)
    {
        this->InitialiseIntVar(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::SetIntVar(const std::string &varName, const int value)
{
    if (m_strIntMap.find(varName) == m_strIntMap.end())
    {
        this->InitialiseIntVar(varName);
        m_strIntMap.at(varName) = value;
        //std::cout<<"TTreeHelper::SetIntVar error m_strIntMap do not have "<<varName<<std::endl;
    }
    else
    {
        
        m_strIntMap.at(varName) = value;
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

int TTreeHelper::GetIntVar(const std::string &varName) const
{
    if (m_strIntMap.find(varName) == m_strIntMap.end())
    {
        std::cout<<"TTreeHelper::GetIntVar error m_strIntMap do not have "<<varName<<std::endl;
    }
    else
    {
        return m_strIntMap.at(varName);
    }
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::InitialiseDoubleVar(const std::string &varName)
{
    if (!m_strDoubleMap.insert(StrIntMap::value_type(varName, 0.f)).second)
        std::cout<<"TTreeHelper::InitialiseDoubleVar error m_strDoubleMap unable to insert "<<varName<<std::endl;
    m_pTTree->Branch(varName.c_str(), &m_strDoubleMap.at(varName), "varName/D");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::InitialiseDoubleVars(const StrSet &varNames)
{
    for (StrSet::const_iterator iter = varNames.begin(), iterEnd = varNames.end(); iter != iterEnd; ++iter)
    {
        this->InitialiseDoubleVar(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TTreeHelper::SetDoubleVar(const std::string &varName, const double value)
{
    if (m_strDoubleMap.find(varName) == m_strDoubleMap.end())
    {
        this->InitialiseDoubleVar(varName);
        m_strDoubleMap.at(varName) = value;
        //std::cout<<"TTreeHelper::SetDoubleVar error m_strDoubleMap do not have "<<varName<<std::endl;
    }
    else
    {
        m_strDoubleMap.at(varName) = value;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

double TTreeHelper::GetDoubleVar(const std::string &varName) const
{
    if (m_strDoubleMap.find(varName) == m_strDoubleMap.end())
    {
        std::cout<<"TTreeHelper::GetDoubleVar error m_strDoubleMap do not have "<<varName<<std::endl;
    }
    else
    {
        return m_strDoubleMap.at(varName);
    }
    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename NumberType> 
NumberType TTreeHelper::GetVar(const std::string &varName) const
{
    if (m_strDoubleMap.find(varName) != m_strDoubleMap.end())
    {
        return m_strDoubleMap.at(varName);
        
    }
    else if (m_strIntMap.find(varName) != m_strIntMap.end())
    {
        return m_strIntMap.at(varName);
    }
    else
    {
        std::cout<<"TTreeHelper::GetVar error m_strDoubleMap or m_strIntMap do not have "<<varName<<std::endl;
    }
    return 0;
}


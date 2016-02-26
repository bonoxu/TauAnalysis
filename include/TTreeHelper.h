#ifndef TTreeHelper_h
#define TTreeHelper_h 1

#include <string>
#include <map>
#include <set>

#include <TFile.h>
#include <TTree.h>

class TTreeHelper 
{

public:
    
    /**
     *  @brief constructor with a tfile and a ttree
     * 
     *  @param tFileName tfile name
     *  @param tFileOption tfile open option
     *  @param tTreeName ttree name
     *  @param tTreeID ttree id
     */
    TTreeHelper(const std::string &tFileName, const std::string &tFileOption, const std::string &tTreeName, const std::string &tTreeID);
    
    /**
     *  @brief destructor
     */
    ~TTreeHelper();
    
    /**
     *  @brief fill node
     */
    void FillNode();
    
    /**
     *  @brief write ttree to tfile
     */
    void WriteTTreeToTFile();
    
    /**
     *  @brief clear double and integer variables
     */
    void Clear();
    
    typedef std::map<std::string, int> StrIntMap;
    typedef std::set<std::string> StrSet;
    
    /**
     *  @brief initialise an integer variable
     * 
     *  @param varName variable name
     */
    void InitialiseIntVar(const std::string &varName);
    
    /**
     *  @brief initialise integer variables
     * 
     *  @param varNames variable names in a set
     */
    void InitialiseIntVars(const StrSet &varNames);
    
    /**
     *  @brief set integer variable
     * 
     *  @param varNames variable name
     *  @param value value
     */
    void SetIntVar(const std::string &varName, const int value);
    
    /**
     *  @brief get integer variable value
     * 
     *  @param varNames variable name
     * 
     *  @return value 
     */
    int GetIntVar(const std::string &varName) const;
    
    typedef std::map<std::string, double> StrDoubleMap;
    
    /**
     *  @brief initialise an double variable
     * 
     *  @param varName variable name
     */
    void InitialiseDoubleVar(const std::string &varName);
    
    /**
     *  @brief initialise double variables
     * 
     *  @param varNames variable names in a set
     */
    void InitialiseDoubleVars(const StrSet &varNames);
    
    /**
     *  @brief set double variable
     * 
     *  @param varNames variable name
     *  @param value value
     */
    void SetDoubleVar(const std::string &varName, const double value);
    
    /**
     *  @brief get double variable value
     * 
     *  @param varNames variable name
     * 
     *  @return value 
     */
    double GetDoubleVar(const std::string &varName) const;
    
protected:
    TTree           *m_pTTree;          // TTree
    TFile           *m_pTFile;          // TFile
    StrIntMap       m_strIntMap;        // string to integer map
    StrDoubleMap    m_strDoubleMap;     // string to double map
};

#endif




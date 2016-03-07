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

    typedef std::map<std::string, int> StrIntMap;
    typedef std::set<std::string> StrSet;
    typedef std::map<std::string, double> StrDoubleMap;
    
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
    
    /**
     *  @brief  set varible
     *
     *  @para	number number
     */
    template <typename NumberType> 
    void SetVar(const std::string &varName, const NumberType number);
    
    /**
     *  @brief  get varible
     *
     *  @param varNames variable name
     *  @para	varName variable name
     */
    template <typename NumberType> 
    NumberType GetVar(const std::string &varName) const;
    
protected:

    /**
     *  @brief  check if the variable is an interger
     *
     *  @para	number number
     * 
     *  @return true if the variable is an interger
     */
    template <typename NumberType> 
    bool IsInt(const NumberType number);

    /**
     *  @brief  check if the variable is a float or double
     *
     *  @para	number number
     * 
     *  @return true if the variable is a float or double
     */
    template <typename NumberType> 
    bool IsFloat(const NumberType number);
    
    TTree           *m_pTTree;          // TTree
    TFile           *m_pTFile;          // TFile
    StrIntMap       m_strIntMap;        // string to integer map
    StrDoubleMap    m_strDoubleMap;     // string to double map
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename NumberType> 
inline bool TTreeHelper::IsInt(const NumberType number)
{
    return (typeid(int) == typeid(number));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename NumberType> 
inline bool TTreeHelper::IsFloat(const NumberType number)
{
    return (typeid(float) == typeid(number) || typeid(double) == typeid(number));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename NumberType> 
inline void TTreeHelper::SetVar(const std::string &varName, const NumberType number)
{
    // ATTN treat all non float/double as int
    return (this->IsFloat(number) ? this->SetDoubleVar(varName, number) : this->SetIntVar(varName, number));
}

#endif




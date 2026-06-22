from json_to_gro_generator.gro_file_compiler import extract_signal_definitions

def prepare_parameters_and_data(sbol_data, parameters):
  """
  Prepares signal definitions, protein actions, and degradation actions
  based on SBOL data and user-defined parameters.

  Args:
    sbol_data (dict): A dictionary containing SBOL data, expected to have keys like
                      "ED", "interactions", "hierarchy", and "components".
    parameters (dict): A dictionary containing user-defined parameters, particularly
                       "signal_parameters" for signal extraction.

  Returns:
    tuple: A tuple containing:
      - signal_definitions (list): List of strings defining signals for GRO.
      - protein_actions_map (dict): Dictionary mapping protein IDs to colors for paint actions.
  """
  signal_definitions = extract_signal_definitions(
    sbol_data["ED"],
    parameters["signal_parameters"]
  )
  
  protein_actions_map = {
    # Green
    "BBa_E0040": "green",   
    "BBa_K2148009": "green", 
    "BBa_K2560042" : "green", 
    "BBa_K2009820": "green", 
    "BBa_K4946001": "green", 
    "BBa_K4159005:": "green",
    "GFP": "green",             
    "gfp": "green",             
    "EGFP": "green",            
    "eGFP": "green",            
    "eGFPuv": "green",          
    "mNeonGreen": "green",      
    "mEmerald": "green",        
    "mGFP": "green",            
    "esmGFP": "green",          
    "mNeonGreen": "green",      
    "mGFP2": "green",           
    "avGFP": "green",

    # Red
    "BBa_E1010": "red",               
    "BBa_K1323009": "red",            
    "BBa_K1688019": "red",             
    "BBa_K3128008": "red",              
    "BBa_K1399000": "red",              
    "BBa_K1399001": "red",               
    "BBa_K1399002": "red",               
    "BBa_K3841014": "red",               
    "RFP": "red",              
    "rfp": "red",               
    "DsRED": "red",             
    "dsRed": "red",             
    "mCherry": "red",           
    "mStrawberry": "red",       
    "mOrange": "orange",        
    "mCherry2": "red",          
    "mRFP1": "red",             
    "mRFP2": "red",             
    "mNeptune": "red",          
    "mRuby": "red",             
    "mKate": "red",             
    "mNeptune2": "red",         
    "mStrawberry2": "red",      
    "mStrawberry3": "red",      
    "mKate2": "red",            
    "mRuby2": "red",        

    # Yellow
    "BBa_E0030": "yellow",  
    "BBa_K592101": "yellow", 
    "BBa_K2656020": "yellow", 
    "BBa_K3427000" : "yellow",
    "BBa_K2656021" : "yellow", 
    "BBa_K592101": "yellow", 
    "BBa_K1323010": "yellow", 
    "BBa_K165005": "yellow", 
    "YFP": "yellow",            
    "yfp": "yellow",           
    "EYFP": "yellow",           
    "eYFP": "yellow",           
    "mCitrine": "yellow",       
    "mVenus": "yellow",         
    "mYFP": "yellow",           
    "mVenus": "yellow",         
    "mVenus2": "yellow",       
    "YFP_LVA" : "yellow",
    "YFP_LOV": "yellow",

    # Cyan 
    "BBa_E0020": "cyan",
    "CFP": "cyan",              
    "cfp": "cyan",              
    "ECFP": "cyan",             
    "eCFP": "cyan",             
    "mTurquoise": "cyan",       
    "mTFP1": "cyan",            
    "mTurq2": "cyan",           
    "mTFP2": "cyan",            
    "mCFP": "cyan",             
    "mTurquoise3": "cyan",      
    "mTFP3": "cyan",           
  }

  return (
    signal_definitions,
    protein_actions_map
  )

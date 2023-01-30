class Utils {

    public static final String ANSI_RESET = "\u001B[0m";
    public static final String ANSI_RED = "\u001B[31m";
    public static final String ANSI_GREEN = "\u001B[32m";
    public static final String ANSI_YELLOW = "\u001B[33m";
    
    static void info(String message, Boolean verbose) {
        if (verbose) {
            println(ANSI_GREEN + "[INFO]\t${message}" + ANSI_RESET)
        }
    }
    static void error(String message) {
        println(ANSI_RED + "[ERROR]\t${message}"+ ANSI_RESET)
        System.exit(1)
    }

    static void warn(String message, Boolean verbose) {
        if (verbose) {
            println(ANSI_YELLOW + "[WARN]\t${message}" + ANSI_RESET)
        }
    }

    static String   dumpParams(filepath,params)    {
        File dumpFile = new File(filepath)
        def date = new Date()
        dumpFile.write "\n//Run parameters ${date}\n"
        dumpFile.append "params {\n"
        params.keySet().sort().each { k -> 
            if (params[k] instanceof String){
                dumpFile.append "${k}='${params[k]}'\n" 
            }
            else {
                dumpFile.append "${k}=${params[k]}\n" 
            }
        }
        dumpFile.append "}\n"
    }

    static void assertFileExists(paramname,path) {
        File file = new File(path)
        if (! file.exists())    {
            this.error("File ${path} does not exists - check your configuration file for paramerter ${paramname}")
        }
    }

}

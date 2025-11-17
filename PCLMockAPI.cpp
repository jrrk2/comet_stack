/**
 * PCLMockAPI.cpp - Implementation of the PCL mock API
 */

#include "PCLMockAPI.h"
#include "PCLThreadMock.h"
#include <pcl/api/APIInterface.h>

#include <map>
#include <mutex>
#include <string>
#include <iostream>
#include <fstream>

namespace pcl_mock {

// Global module handle
static void* g_module_handle = nullptr;

// Function mapping
static std::map<std::string, void*> g_function_map;
static std::mutex g_function_map_mutex;

// Logging settings
static bool g_debug_logging = false;
static std::ofstream g_log_file;

// Log a debug message
void LogDebug(const std::string& message) {
    if (!g_debug_logging) return;
    
    if (g_log_file.is_open()) {
        g_log_file << message << std::endl;
    } else {
        std::cout << "[PCLMockAPI] " << message << std::endl;
    }
}

    // Create a stub function for missing functions
    // We'll use a simple global function that just returns nullptr
    static void* unimplemented_function(void) {
        LogDebug("Called unimplemented function");
        return nullptr;
    }
      
// Function resolver implementation
void* mock_function_resolver(const char* name) {
    if (!name) return nullptr;
    
    std::string func_name = name;
    LogDebug("Resolving function: " + func_name);
    
    std::lock_guard<std::mutex> lock(g_function_map_mutex);
    auto it = g_function_map.find(func_name);
    if (it != g_function_map.end()) {
        LogDebug("Found implementation for: " + func_name);
        return it->second;
    }
    
    // Create a default handler for missing functions
    LogDebug("No implementation found for: " + func_name);
    
    // Store this handler so we don't create a new one each time
    g_function_map[func_name] = (void*)unimplemented_function;
    
    return g_function_map[func_name];
}

// Register a function with the mock API
void RegisterFunction(const char* name, void* func) {
    if (!name || !func) return;
    
    std::lock_guard<std::mutex> lock(g_function_map_mutex);
    g_function_map[name] = func;
    LogDebug("Registered function: " + std::string(name));
}

// Initialize the mock API
void InitializeMockAPI() {
    // Register thread-related functions
    RegisterThreadFunctions();
    RegisterGlobalFunctions();
    RegisterUIFunctions();

    // Create API interface with mock function resolver if not already created
    if (!API) {
      API = new pcl::APIInterface(pcl_mock::GetMockFunctionResolver());
    }

    // Set module handle if not already set
    if (!GetModuleHandle()) {
      pcl_mock::SetModuleHandle((void*)0x12345678);
    }
    
    LogDebug("Mock API initialized");
}

// Get the function resolver
function_resolver GetMockFunctionResolver() {
    return mock_function_resolver;
}

// Set the module handle
void SetModuleHandle(void* handle) {
    g_module_handle = handle;
    LogDebug("Module handle set to: " + std::to_string((uintptr_t)handle));
}

// Get the module handle
void* GetModuleHandle() {
    return g_module_handle;
}

// Enable or disable debug logging
void SetDebugLogging(bool enabled) {
    g_debug_logging = enabled;
}

// Set the log file
void SetLogFile(const std::string& filename) {
    if (g_log_file.is_open()) {
        g_log_file.close();
    }
    
    if (!filename.empty()) {
        g_log_file.open(filename);
        if (!g_log_file.is_open()) {
            std::cerr << "Failed to open log file: " << filename << std::endl;
        }
    }
}

// Mock implementations for Global API functions

// Simple error code storage
static int g_last_error = 0;

// Mock console handle
static void* g_console_handle = (void*)0xDEADBEEF;

// Mock pixel traits LUT - create a simple lookup table
struct MockPixelTraitsLUT {
    // These would normally be function pointers, but we'll make them simple values
    int sample_format;
    int bytes_per_sample;
    int bits_per_sample;
    double min_sample_value;
    double max_sample_value;
};

// Create mock LUTs for different pixel formats
static MockPixelTraitsLUT g_pixel_luts[16] = {
    // Format 0: 8-bit unsigned integer
    { 0, 1, 8, 0.0, 255.0 },
    // Format 1: 16-bit unsigned integer  
    { 1, 2, 16, 0.0, 65535.0 },
    // Format 2: 32-bit unsigned integer
    { 2, 4, 32, 0.0, 4294967295.0 },
    // Format 3: 32-bit IEEE 754 floating point
    { 3, 4, 32, 0.0, 1.0 },
    // Format 4: 64-bit IEEE 754 floating point
    { 4, 8, 64, 0.0, 1.0 },
    // Add more formats as needed...
};

// Mock for GetPixelTraitsLUT
void* GetPixelTraitsLUT(int format) {
    LogDebug("GetPixelTraitsLUT called with format: " + std::to_string(format));
    
    if (format >= 0 && format < 16) {
        return &g_pixel_luts[format];
    }
    return &g_pixel_luts[0]; // Default to format 0
}

// Mock for GetConsole
void* GetConsole() {
    LogDebug("GetConsole called");
    return g_console_handle;
}

// Mock for LastError
int LastError() {
    LogDebug("LastError called, returning: " + std::to_string(g_last_error));
    return g_last_error;
}

// Mock for setting error
void SetLastError(int error_code) {
    g_last_error = error_code;
    LogDebug("SetLastError called with: " + std::to_string(error_code));
}

// Mock for ClearError
void ClearError() {
    g_last_error = 0;
    LogDebug("ClearError called");
}

// Mock for ProcessEvents
int ProcessEvents() {
    LogDebug("ProcessEvents called");
    return 1; // Success
}

// Mock for GetApplicationInstanceSlot
int GetApplicationInstanceSlot() {
    LogDebug("GetApplicationInstanceSlot called");
    return 0; // Root slot
}

// Mock for GetProcessStatus
uint32_t GetProcessStatus() {
    LogDebug("GetProcessStatus called");
    // Return a status that indicates not aborted (bit 31 clear)
    // PCL checks if bit 31 (0x80000000) is set to determine if process should abort
    return 0x00000000; // Normal status, not aborted
}

// Mock for WriteConsole
int WriteConsole(void* console_handle, const char* text, int append_newline) {
    LogDebug("WriteConsole called with text: " + std::string(text ? text : "(null)"));
    
    if (!text) {
        return 0; // api_false
    }
    
    // Just output to stdout for now
    std::cout << text;
    if (append_newline) {
        std::cout << std::endl;
    }
    
    return 1; // api_true - success
}


// Mock for GetGlobalFlag
int GetGlobalFlag(const char* flag_name, int* value) {
    LogDebug("GetGlobalFlag called with: " + std::string(flag_name ? flag_name : "(null)"));
    
    if (!value) {
        return 0; // api_false
    }
    
    // Return default values for common flags
    if (flag_name) {
        std::string name = flag_name;
        if (name.find("Abort") != std::string::npos) {
            *value = 0; // Not aborted
        } else if (name.find("Debug") != std::string::npos) {
            *value = 0; // Debug off
        } else {
            *value = 0; // Default to false/off
        }
    } else {
        *value = 0;
    }
    
    return 1; // api_true - success
}

// Mock for GetGlobalInteger
int GetGlobalInteger(const char* int_name, int* value, int create_if_not_exists) {
    LogDebug("GetGlobalInteger called with: " + std::string(int_name ? int_name : "(null)"));
    
    if (!value) {
        return 0; // api_false
    }
    
    // Return default values for common integers
    if (int_name) {
        std::string name = int_name;
        if (name.find("MaxProcessors") != std::string::npos) {
            *value = 4; // Default to 4 processors
        } else if (name.find("ThreadPriority") != std::string::npos) {
            *value = 3; // Normal priority
        } else {
            *value = 0; // Default value
        }
    } else {
        *value = 0;
    }
    
    return 1; // api_true - success
}

// Mock for GetUIObjectRefCount
int GetUIObjectRefCount(void* ui_object) {
    LogDebug("GetUIObjectRefCount called with object: " + std::to_string((uintptr_t)ui_object));
    
    if (!ui_object) {
        return 0; // No references for null object
    }
    
    // Return 1 to indicate the object exists and has at least one reference
    return 1;
}

// Mock for DetachFromUIObject
int DetachFromUIObject(void* module_handle, void* ui_object) {
    LogDebug("DetachFromUIObject called with module: " + std::to_string((uintptr_t)module_handle) + 
             ", object: " + std::to_string((uintptr_t)ui_object));
    
    if (!ui_object) {
        return 0; // api_false - can't detach from null object
    }
    
    // Always succeed in detaching (in a real implementation, this would decrease reference count)
    return 1; // api_true - success
}

// Add these to your RegisterGlobalFunctions() - create this new function
void RegisterGlobalFunctions() {
    RegisterFunction("Global/GetPixelTraitsLUT", (void*)GetPixelTraitsLUT);
    RegisterFunction("Global/GetConsole", (void*)GetConsole);
    RegisterFunction("Global/LastError", (void*)LastError);
    RegisterFunction("Global/ClearError", (void*)ClearError);
    
    // Add some other commonly used Global functions
    RegisterFunction("Global/ProcessEvents", (void*)ProcessEvents);
    RegisterFunction("Global/GetApplicationInstanceSlot", (void*)GetApplicationInstanceSlot);
    RegisterFunction("Global/GetProcessStatus", (void*)GetProcessStatus);
    RegisterFunction("Global/WriteConsole", (void*)WriteConsole);
    RegisterFunction("Global/GetGlobalFlag", (void*)GetGlobalFlag);
    RegisterFunction("Global/GetGlobalInteger", (void*)GetGlobalInteger);
}

// Register UI functions
void RegisterUIFunctions() {
    RegisterFunction("UI/GetUIObjectRefCount", (void*)GetUIObjectRefCount);
    RegisterFunction("UI/DetachFromUIObject", (void*)DetachFromUIObject);
}

} // namespace pcl_mock

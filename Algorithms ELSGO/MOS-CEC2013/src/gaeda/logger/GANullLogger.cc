#include "GANullLogger.h"

GANullLogger::GANullLogger (char* path, unsigned int stats_display_freq, GALogger::LogLevel log_level) : GALogger (stats_display_freq, log_level) {}
GANullLogger::~GANullLogger() {}

/* ----------------------------------------------------------------------------
   libconfig - A structured configuration file parsing library
   Copyright (C) 2005-2007  Mark A Lindner
 
   This file is part of libconfig.
    
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   as published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.
    
   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
    
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
   ----------------------------------------------------------------------------
*/

#include "libconfigcpp.h"

#ifdef _MSC_VER
#pragma warning (disable: 4996)
#endif

using namespace libconfig;

// ---------------------------------------------------------------------------

static int __toTypeCode(Setting::Type type)
{
  int typecode;

  switch(type)
  {
    case Setting::TypeGroup:
      typecode = CONFIG_TYPE_GROUP;
      break;
      
    case Setting::TypeInt:
      typecode = CONFIG_TYPE_INT;
      break;

    case Setting::TypeFloat:
      typecode = CONFIG_TYPE_FLOAT;
      break;

    case Setting::TypeString:
      typecode = CONFIG_TYPE_STRING;
      break;

    case Setting::TypeBoolean:
      typecode = CONFIG_TYPE_BOOL;
      break;

    case Setting::TypeArray:
      typecode = CONFIG_TYPE_ARRAY;
      break;

    case Setting::TypeList:
      typecode = CONFIG_TYPE_LIST;
      break;
      
    default:
      typecode = CONFIG_TYPE_NONE;
  }

  return(typecode);
}

// ---------------------------------------------------------------------------

void Config::ConfigDestructor(void *arg)
{
  delete reinterpret_cast<Setting *>(arg);
}

// ---------------------------------------------------------------------------

Config::Config()
{
  config_init(& _config);
  config_set_destructor(& _config, ConfigDestructor);
}

// ---------------------------------------------------------------------------

Config::~Config()
{
  config_destroy(& _config);
}

// ---------------------------------------------------------------------------

void Config::setAutoConvert(bool flag)
{
  config_set_auto_convert(& _config, (flag ? CONFIG_TRUE : CONFIG_FALSE));
}

// ---------------------------------------------------------------------------

bool Config::getAutoConvert() const
{
  return(config_get_auto_convert(& _config) != CONFIG_FALSE);
}

// ---------------------------------------------------------------------------

void Config::read(FILE *stream) throw(ParseException)
{
  if(! config_read(& _config, stream))
    throw ParseException(config_error_line(& _config),
                         config_error_text(& _config));
}

// ---------------------------------------------------------------------------

void Config::write(FILE *stream) const
{
  config_write(& _config, stream);
}

// ---------------------------------------------------------------------------

void Config::readFile(const char *filename) throw(FileIOException,
                                                  ParseException)
{
  FILE *f = fopen(filename, "rt");
  if(f == NULL)
    throw FileIOException();
  try
  {
    read(f);
    fclose(f);
  }
  catch(ParseException& p)
  {
    fclose(f);
    throw p;
  }
}

// ---------------------------------------------------------------------------

void Config::writeFile(const char *filename) throw(FileIOException)
{
  if(! config_write_file(& _config, filename))
    throw FileIOException();
}

// ---------------------------------------------------------------------------

Setting & Config::lookup(const char *path) const
  throw(SettingNotFoundException)
{
  config_setting_t *s = config_lookup(& _config, path);
  if(! s)
    throw SettingNotFoundException();

  return(Setting::wrapSetting(s));
}

// ---------------------------------------------------------------------------

bool Config::exists(const char *path) const throw()
{
  config_setting_t *s = config_lookup(& _config, path);

  return(s != NULL);
}

// ---------------------------------------------------------------------------

#define CONFIG_LOOKUP_NO_EXCEPTIONS(P, T, V)    \
  try                                           \
  {                                             \
    Setting &s = lookup(P);                     \
    V = (T)s;                                   \
    return(true);                               \
  }                                             \
  catch(ConfigException)                        \
  {                                             \
    return(false);                              \
  }

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, bool &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, bool, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, long &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, long, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, unsigned long &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, unsigned long, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, int &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, int, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, unsigned int &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, unsigned int, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, double &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, double, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, float &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, float, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, const char *&value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, const char *, value);
}

// ---------------------------------------------------------------------------

bool Config::lookupValue(const char *path, std::string &value) const throw()
{
  CONFIG_LOOKUP_NO_EXCEPTIONS(path, const char *, value);
}

// ---------------------------------------------------------------------------

Setting & Config::getRoot() const
{
  return(Setting::wrapSetting(config_root_setting(& _config)));
}

// ---------------------------------------------------------------------------

Setting::Setting(config_setting_t *setting)
  : _setting(setting)
{
  switch(config_setting_type(setting))
  {
    case CONFIG_TYPE_GROUP:
      _type = TypeGroup;
      break;

    case CONFIG_TYPE_INT:
      _type = TypeInt;
      break;

    case CONFIG_TYPE_FLOAT:
      _type = TypeFloat;
      break;

    case CONFIG_TYPE_STRING:
      _type = TypeString;
      break;

    case CONFIG_TYPE_BOOL:
      _type = TypeBoolean;
      break;

    case CONFIG_TYPE_ARRAY:
      _type = TypeArray;
      break;

    case CONFIG_TYPE_LIST:
      _type = TypeList;
      break;

    case CONFIG_TYPE_NONE:
    default:
      _type = TypeNone;
      break;
  }

  switch(config_setting_get_format(setting))
  {
    case CONFIG_FORMAT_HEX:
      _format = FormatHex;
      break;

    case CONFIG_FORMAT_DEFAULT:
    default:
      _format = FormatDefault;
      break;
  }
}

// ---------------------------------------------------------------------------

Setting::~Setting() throw()
{
  _setting = NULL;
}

// ---------------------------------------------------------------------------

void Setting::setFormat(Format format) throw()
{
  if(_type == TypeInt)
  {
    if(format ==FormatHex)
      _format = FormatHex;
    else
      _format = FormatDefault;
  }
  else
    _format = FormatDefault;
}

// ---------------------------------------------------------------------------

Setting::operator bool() const throw(SettingTypeException) 
{
  assertType(TypeBoolean);

  return(config_setting_get_bool(_setting) ? true : false);
}

// ---------------------------------------------------------------------------

Setting::operator long() const throw(SettingTypeException)
{
  assertType(TypeInt);

  return(config_setting_get_int(_setting));
}

// ---------------------------------------------------------------------------

Setting::operator unsigned long() const throw(SettingTypeException)
{
  assertType(TypeInt);

  long v = config_setting_get_int(_setting);

  if(v < 0)
    v = 0;
  
  return(static_cast<unsigned long>(v));
}

// ---------------------------------------------------------------------------

Setting::operator int() const throw(SettingTypeException)
{
  assertType(TypeInt);

  // may cause loss of precision:
  return(static_cast<int>(config_setting_get_int(_setting)));
}

// ---------------------------------------------------------------------------

Setting::operator unsigned int() const throw(SettingTypeException)
{
  assertType(TypeInt);

  long v = config_setting_get_int(_setting);

  if(v < 0)
    v = 0;

  return(static_cast<unsigned int>(v));
}

// ---------------------------------------------------------------------------

Setting::operator double() const throw(SettingTypeException)
{
  assertType(TypeFloat);

  return(config_setting_get_float(_setting));
}

// ---------------------------------------------------------------------------

Setting::operator float() const throw(SettingTypeException)
{
  assertType(TypeFloat);

  // may cause loss of precision:
  return(static_cast<float>(config_setting_get_float(_setting)));
}

// ---------------------------------------------------------------------------

Setting::operator const char *() const throw(SettingTypeException)
{
  assertType(TypeString);

  return(config_setting_get_string(_setting));
}

// ---------------------------------------------------------------------------

Setting::operator std::string() const throw(SettingTypeException)
{
  assertType(TypeString);

  const char *s = config_setting_get_string(_setting);

  std::string str;
  if(s)
    str = s;

  return(str);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator=(bool value) throw(SettingTypeException)
{
  assertType(TypeBoolean);

  config_setting_set_bool(_setting, value);

  return(*this);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator=(long value) throw(SettingTypeException)
{
  assertType(TypeInt);

  config_setting_set_int(_setting, value);

  return(*this);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator=(int value) throw(SettingTypeException)
{
  assertType(TypeInt);

  long cvalue = static_cast<long>(value);
  
  config_setting_set_int(_setting, cvalue);

  return(*this);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator=(const double &value) throw(SettingTypeException)
{
  assertType(TypeFloat);

  config_setting_set_float(_setting, value);

  return(*this);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator=(float value) throw(SettingTypeException)
{
  assertType(TypeFloat);

  double cvalue = static_cast<double>(value);

  config_setting_set_float(_setting, cvalue);

  return(*this);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator=(const char *value) throw(SettingTypeException)
{
  assertType(TypeString);

  config_setting_set_string(_setting, value);

  return(*this);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator=(const std::string &value)
  throw(SettingTypeException)
{
  assertType(TypeString);

  config_setting_set_string(_setting, value.c_str());

  return(*this);
}

// ---------------------------------------------------------------------------

Setting & Setting::operator[](int i) const
  throw(SettingTypeException, SettingNotFoundException)
{
  if((_type != TypeArray) && (_type != TypeGroup) && (_type != TypeList))
    throw SettingTypeException();
  
  config_setting_t *setting = config_setting_get_elem(_setting, i);

  if(! setting)
    throw SettingNotFoundException();

  return(wrapSetting(setting));
}

// ---------------------------------------------------------------------------

Setting & Setting::operator[](const char *key) const
  throw(SettingTypeException, SettingNotFoundException)
{
  assertType(TypeGroup);

  config_setting_t *setting = config_setting_get_member(_setting, key);

  if(! setting)
    throw SettingNotFoundException();

  return(wrapSetting(setting));
}

// ---------------------------------------------------------------------------

#define SETTING_LOOKUP_NO_EXCEPTIONS(K, T, V)   \
  try                                           \
  {                                             \
    Setting &s = operator[](K);                 \
    V = (T)s;                                   \
    return(true);                               \
  }                                             \
  catch(ConfigException)                        \
  {                                             \
    return(false);                              \
  }

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, bool &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, bool, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, long &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, long, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, unsigned long &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, unsigned long, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, int &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, int, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, unsigned int &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, unsigned int, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, double &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, double, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, float &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, float, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, const char *&value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, const char *, value);
}

// ---------------------------------------------------------------------------

bool Setting::lookupValue(const char *name, std::string &value) const throw()
{
  SETTING_LOOKUP_NO_EXCEPTIONS(name, const char *, value);
}

// ---------------------------------------------------------------------------

bool Setting::exists(const char *name) const throw()
{
  if(_type != TypeGroup)
    return(false);

  config_setting_t *setting = config_setting_get_member(_setting, name);

  return(setting != NULL);
}

// ---------------------------------------------------------------------------

int Setting::getLength() const throw()
{
  return(config_setting_length(_setting));
}

// ---------------------------------------------------------------------------

const char * Setting::getName() const throw()
{
  return(config_setting_name(_setting));
}

// ---------------------------------------------------------------------------

void Setting::remove(const char *name)
  throw(SettingTypeException, SettingNotFoundException)
{
  assertType(TypeGroup);
  
  if(! config_setting_remove(_setting, name))
    throw SettingNotFoundException();
}

// ---------------------------------------------------------------------------

Setting & Setting::add(const char *name, Setting::Type type)
  throw(SettingTypeException, SettingExistsException)
{
  assertType(TypeGroup);
  
  int typecode = __toTypeCode(type);

  if(typecode == CONFIG_TYPE_NONE)
    throw SettingTypeException();

  config_setting_t *setting = config_setting_add(_setting, name, typecode);

  if(! setting)
    throw SettingExistsException();

  return(wrapSetting(setting));
}

// ---------------------------------------------------------------------------

Setting & Setting::add(Setting::Type type) throw(SettingTypeException)
{
  if((_type != TypeArray) && (_type != TypeList))
    throw SettingTypeException();

  if(_type == TypeArray)
  {  
    if(getLength() > 0)
    {
      Setting::Type atype = operator[](0).getType();
      if(type != atype)
        throw SettingTypeException();
    }
    else
    {
      if((type != TypeInt) && (type != TypeFloat) && (type != TypeString)
         && (type != TypeBoolean))
        throw SettingTypeException();
    }
  }

  int typecode = __toTypeCode(type);
  config_setting_t *s = config_setting_add(_setting, NULL, typecode);

  Setting &ns = wrapSetting(s);

  switch(type)
  {
    case TypeInt:
      ns = 0;
      break;

    case TypeFloat:
      ns = 0.0;
      break;

    case TypeString:
      ns = (char *)NULL;
      break;

    case TypeBoolean:
      ns = false;
      break;

    default:
      // won't happen
      break;
  }

  return(ns);
}

// ---------------------------------------------------------------------------

void Setting::assertType(Setting::Type type) const throw(SettingTypeException)
{
  if(type != _type)
  {
    if(!(isNumber() && config_get_auto_convert(_setting->config)
         && ((type == TypeInt) || (type == TypeFloat))))
      throw SettingTypeException();
  }
}

// ---------------------------------------------------------------------------

Setting & Setting::wrapSetting(config_setting_t *s)
{
  Setting *setting = NULL;
  
  void *hook = config_setting_get_hook(s);
  if(! hook)
  {
    setting = new Setting(s);
    config_setting_set_hook(s, reinterpret_cast<void *>(setting));
  }
  else
    setting = reinterpret_cast<Setting *>(hook);

  return(*setting);
}

// ---------------------------------------------------------------------------
// eof

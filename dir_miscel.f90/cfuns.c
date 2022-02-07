#if defined(OSF1)
int putenv_(char *str)
{
  return putenv(str);
}
#else
void dumdum(void)
{
;
}
#endif

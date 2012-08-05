/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef QUIP_CRC64
#define QUIP_CRC64

#include <stdlib.h>
#include <stdint.h>

uint64_t crc64_update(uint8_t* buf, size_t bufsize, uint64_t crc);

#endif

/*
 * Copyright (c) 2010, University of Twente
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in
 *   the documentation and/or other materials provided with the
 *   distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*
 *  speedklein80.h
 *
 *  Description: klein-80 with a clear structure.
 *  Created on: 2010-4-7
 *  Last modified: 2010-8-18
 *  Author: Zheng Gong
 */

#include "kleinSbox.h"
#ifndef KLEIN80_H_
#define KLEIN80_H_

#ifdef __cplusplus
extern "C" {
#endif

//rounds of Klein-80;
#define ROUNDS_80 16

#define klein80_encrypt(plain, key, cipher) klein80_encrypt_rounds((plain), (key), ROUNDS_80, (cipher))
#define klein80_decrypt(cipher, key, plain) klein80_decrypt_rounds((cipher), (key), ROUNDS_80, (plain))

void klein80_encrypt_rounds(const uint8_t *plain, const uint8_t *key, const uint8_t rounds, uint8_t *cipher);

void klein80_decrypt_rounds(const uint8_t *cipher, const uint8_t *key, const uint8_t rounds, uint8_t *plain);

#ifdef __cplusplus
}
#endif

#endif /* KLEIN80_H_ */

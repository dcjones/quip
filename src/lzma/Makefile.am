
noinst_LTLIBRARIES = liblzma.la

liblzma_la_SOURCES = \
	mythread.h sysdefs.h \
	tuklib_common.h tuklib_config.h tuklib_integer.h \
	tuklib_physmem.h tuklib_physmem.c

liblzma_la_SOURCES += \
	lzma/lzma_common.h \
	lzma/fastpos.h \
	lzma/fastpos_table.c \
	lzma/lzma_encoder.h \
	lzma/lzma_encoder.c \
	lzma/lzma_encoder_presets.c \
	lzma/lzma_encoder_private.h \
	lzma/lzma_encoder_optimum_fast.c \
	lzma/lzma_encoder_optimum_normal.c \
	lzma/lzma_decoder.h \
	lzma/lzma_decoder.c \
	lzma/lzma2_encoder.h \
	lzma/lzma2_encoder.c \
	lzma/lzma2_decoder.h \
	lzma/lzma2_decoder.c

liblzma_la_SOURCES += \
	lz/lz_decoder.h \
	lz/lz_decoder.c \
	lz/lz_encoder.h \
	lz/lz_encoder.c \
	lz/lz_encoder_hash.h \
	lz/lz_encoder_hash_table.h \
	lz/lz_encoder_mf.c

liblzma_la_SOURCES += \
	rangecoder/range_encoder.h \
	rangecoder/price.h \
	rangecoder/price_table.c \
	rangecoder/range_decoder.h

liblzma_la_SOURCES += \
	check/check.h \
	check/check.c \
	check/crc_macros.h \
	check/crc32_table.c \
	check/crc32_table_le.h \
	check/crc32_table_be.h \
	check/crc32_fast.c \
	check/crc64_table.c \
	check/crc64_table_le.h \
	check/crc64_table_be.h \
	check/crc64_fast.c \
	check/sha256.c

liblzma_la_SOURCES += \
	simple/simple_coder.h \
	simple/simple_decoder.h \
	simple/simple_encoder.h

liblzma_la_SOURCES += \
	delta/delta_common.h \
	delta/delta_decoder.h \
	delta/delta_encoder.h

liblzma_la_SOURCES += \
	common/common.c \
	common/common.h \
	common/block_util.c \
	common/easy_preset.c \
	common/easy_preset.h \
	common/filter_common.c \
	common/filter_common.h \
	common/hardware_physmem.c \
	common/index.c \
	common/index.h \
	common/stream_flags_common.c \
	common/stream_flags_common.h \
	common/vli_size.c \
	common/alone_encoder.c \
	common/block_buffer_encoder.c \
	common/block_encoder.c \
	common/block_encoder.h \
	common/block_header_encoder.c \
	common/easy_buffer_encoder.c \
	common/easy_encoder.c \
	common/easy_encoder_memusage.c \
	common/filter_buffer_encoder.c \
	common/filter_encoder.c \
	common/filter_encoder.h \
	common/filter_flags_encoder.c \
	common/index_encoder.c \
	common/index_encoder.h \
	common/stream_buffer_encoder.c \
	common/stream_encoder.c \
	common/stream_encoder.h \
	common/stream_flags_encoder.c \
	common/vli_encoder.c \
	common/alone_decoder.c \
	common/alone_decoder.h \
	common/auto_decoder.c \
	common/block_buffer_decoder.c \
	common/block_decoder.c \
	common/block_decoder.h \
	common/block_header_decoder.c \
	common/easy_decoder_memusage.c \
	common/filter_buffer_decoder.c \
	common/filter_decoder.c \
	common/filter_decoder.h \
	common/filter_flags_decoder.c \
	common/index_decoder.c \
	common/index_hash.c \
	common/stream_buffer_decoder.c \
	common/stream_decoder.c \
	common/stream_decoder.h \
	common/stream_flags_decoder.c \
	common/vli_decoder.c

liblzma_la_SOURCES += \
	api/lzma.h \
	api/lzma/base.h \
	api/lzma/bcj.h \
	api/lzma/block.h \
	api/lzma/check.h \
	api/lzma/container.h \
	api/lzma/delta.h \
	api/lzma/filter.h \
	api/lzma/hardware.h \
	api/lzma/index.h \
	api/lzma/index_hash.h \
	api/lzma/lzma.h \
	api/lzma/stream_flags.h \
	api/lzma/version.h \
	api/lzma/vli.h

liblzma_la_CPPFLAGS = -Iapi -Icommon -Ilz -Irangecoder -Icheck -Ilzma -Isimple -Idelta

liblzma_la_CFLAGS = -w



import phangsPipeline as ppl
from phangsPipeline.handlerImagingChunked import ImagingChunkedHandler

# Quick test cases for handlerChunkedImaging, not currently set up as proper tests

master_key_file = "/path/to/master_key.txt"

target = "ngc0253"
interf_config = "12m"
line_product = "co21"

imaging_method = "tclean"  # or sdintimaging

chunksize = 10
chunk_temp_path = "/path/to/temp_dir"

key_handler = ppl.KeyHandler(master_key=master_key_file)

high_snr = 2.0
low_snr = 1.0
absolute = True

convergence_fracflux = 0.001
singlescale_threshold_value = 1

# Case 1) Run all chunks at once, don't use a temp directory

im_handler = ImagingChunkedHandler(
                        target=target,
                        config=interf_config,
                        product=line_product,
                        key_handler=key_handler,
                        imaging_method=imaging_method,
                        chunksize=chunksize,
                        make_temp_dir=False,
                        temp_path=chunk_temp_path,
                        copy_ms_to_temp=True,
                    )

im_handler.run_imaging(chunk_num=None,
                       do_dirty_image=True,
                       do_revert_to_dirty=True,
                       do_read_clean_mask=True,
                       do_multiscale_clean=True,
                       do_revert_to_multiscale=True,
                       do_singlescale_clean=True,
                       do_revert_to_singlescale=True,
                       do_singlescale_mask=True,
                       singlescale_mask_absolute=absolute,
                       singlescale_mask_high_snr=high_snr,
                       singlescale_mask_low_snr=low_snr,
                       convergence_fracflux=convergence_fracflux,
                       singlescale_threshold_value=singlescale_threshold_value,
                       do_recombine_cubes=True,
                       do_export_to_fits=True,
                       do_cleanup=True,
                       overwrite=True,
                       )

# Case 2) Run all chunks at once, use a temp dir, don't copy MS

im_handler = ImagingChunkedHandler(
                        target=target,
                        config=interf_config,
                        product=line_product,
                        key_handler=key_handler,
                        imaging_method=imaging_method,
                        chunksize=chunksize,
                        make_temp_dir=True,
                        temp_path=chunk_temp_path,
                        copy_ms_to_temp=False,
                    )

im_handler.run_imaging(chunk_num=None,
                       do_dirty_image=True,
                       do_revert_to_dirty=True,
                       do_read_clean_mask=True,
                       do_multiscale_clean=True,
                       do_revert_to_multiscale=True,
                       do_singlescale_clean=True,
                       do_revert_to_singlescale=True,
                       do_singlescale_mask=True,
                       singlescale_mask_absolute=absolute,
                       singlescale_mask_high_snr=high_snr,
                       singlescale_mask_low_snr=low_snr,
                       convergence_fracflux=convergence_fracflux,
                       singlescale_threshold_value=singlescale_threshold_value,
                       do_recombine_cubes=True,
                       do_export_to_fits=True,
                       do_cleanup=True,
                       overwrite=True,
                       )

# Case 3) Run all chunks at once, use a temp dir, copy MS

im_handler = ImagingChunkedHandler(
                        target=target,
                        config=interf_config,
                        product=line_product,
                        key_handler=key_handler,
                        imaging_method=imaging_method,
                        chunksize=chunksize,
                        make_temp_dir=True,
                        temp_path=chunk_temp_path,
                        copy_ms_to_temp=True,
                    )

im_handler.run_imaging(chunk_num=None,
                       do_dirty_image=True,
                       do_revert_to_dirty=True,
                       do_read_clean_mask=True,
                       do_multiscale_clean=True,
                       do_revert_to_multiscale=True,
                       do_singlescale_clean=True,
                       do_revert_to_singlescale=True,
                       do_singlescale_mask=True,
                       singlescale_mask_absolute=absolute,
                       singlescale_mask_high_snr=high_snr,
                       singlescale_mask_low_snr=low_snr,
                       convergence_fracflux=convergence_fracflux,
                       singlescale_threshold_value=singlescale_threshold_value,
                       do_recombine_cubes=False,
                       do_export_to_fits=True,
                       do_cleanup=True,
                       overwrite=True,
                       )

# Case 4) Run chunks as a loop (or farm off to other machines), don't use a temp dir

im_handler = ImagingChunkedHandler(
                        target=target,
                        config=interf_config,
                        product=line_product,
                        key_handler=key_handler,
                        imaging_method=imaging_method,
                        chunksize=chunksize,
                        make_temp_dir=False,
                        temp_path=chunk_temp_path,
                        copy_ms_to_temp=True,
                    )

for chunk_num in range(im_handler.nchunks):

    im_handler.run_imaging(chunk_num=chunk_num,
                           do_dirty_image=True,
                           do_revert_to_dirty=True,
                           do_read_clean_mask=True,
                           do_multiscale_clean=True,
                           do_revert_to_multiscale=True,
                           do_singlescale_clean=True,
                           do_revert_to_singlescale=True,
                           do_singlescale_mask=True,
                           singlescale_mask_absolute=absolute,
                           singlescale_mask_high_snr=high_snr,
                           singlescale_mask_low_snr=low_snr,
                           convergence_fracflux=convergence_fracflux,
                           singlescale_threshold_value=singlescale_threshold_value,
                           do_recombine_cubes=False,
                           do_export_to_fits=True,
                           do_cleanup=True,
                           overwrite=True,
                           )

im_handler.task_complete_gather_into_cubes()
im_handler.task_export_to_fits()
im_handler.task_cleanup()

# Case 5) Run chunks as a loop (or farm off to other machines), use temp dir, don't copy MS

im_handler = ImagingChunkedHandler(
    target=target,
    config=interf_config,
    product=line_product,
    key_handler=key_handler,
    imaging_method=imaging_method,
    chunksize=chunksize,
    make_temp_dir=True,
    temp_path=chunk_temp_path,
    copy_ms_to_temp=False,
)

for chunk_num in range(im_handler.nchunks):
    im_handler.run_imaging(chunk_num=chunk_num,
                           do_dirty_image=True,
                           do_revert_to_dirty=True,
                           do_read_clean_mask=True,
                           do_multiscale_clean=True,
                           do_revert_to_multiscale=True,
                           do_singlescale_clean=True,
                           do_revert_to_singlescale=True,
                           do_singlescale_mask=True,
                           singlescale_mask_absolute=absolute,
                           singlescale_mask_high_snr=high_snr,
                           singlescale_mask_low_snr=low_snr,
                           convergence_fracflux=convergence_fracflux,
                           singlescale_threshold_value=singlescale_threshold_value,
                           do_recombine_cubes=False,
                           do_export_to_fits=True,
                           do_cleanup=True,
                           overwrite=True,
                           )

im_handler.task_complete_gather_into_cubes()
im_handler.task_export_to_fits()
im_handler.task_cleanup()

# Case 6) Run chunks as a loop (or farm off to other machines), use a temp dir, copy MS

im_handler = ImagingChunkedHandler(
    target=target,
    config=interf_config,
    product=line_product,
    key_handler=key_handler,
    imaging_method=imaging_method,
    chunksize=chunksize,
    make_temp_dir=True,
    temp_path=chunk_temp_path,
    copy_ms_to_temp=True,
)

for chunk_num in range(im_handler.nchunks):
    im_handler.run_imaging(chunk_num=chunk_num,
                           do_dirty_image=True,
                           do_revert_to_dirty=True,
                           do_read_clean_mask=True,
                           do_multiscale_clean=True,
                           do_revert_to_multiscale=True,
                           do_singlescale_clean=True,
                           do_revert_to_singlescale=True,
                           do_singlescale_mask=True,
                           singlescale_mask_absolute=absolute,
                           singlescale_mask_high_snr=high_snr,
                           singlescale_mask_low_snr=low_snr,
                           convergence_fracflux=convergence_fracflux,
                           singlescale_threshold_value=singlescale_threshold_value,
                           do_recombine_cubes=True,
                           do_export_to_fits=True,
                           do_cleanup=True,
                           overwrite=True,
                           )

im_handler.task_complete_gather_into_cubes()
im_handler.task_export_to_fits()
im_handler.task_cleanup()

# Case 7) Run chunks as a loop, use a temp dir, mimic as on a cluster so re-instantiate the handler after

im_handler = ImagingChunkedHandler(
    target=target,
    config=interf_config,
    product=line_product,
    key_handler=key_handler,
    imaging_method=imaging_method,
    chunksize=chunksize,
    make_temp_dir=True,
    temp_path=chunk_temp_path,
    copy_ms_to_temp=False,
)

for chunk_num in range(im_handler.nchunks):
    im_handler.run_imaging(chunk_num=chunk_num,
                           do_dirty_image=True,
                           do_revert_to_dirty=True,
                           do_read_clean_mask=True,
                           do_multiscale_clean=True,
                           do_revert_to_multiscale=True,
                           do_singlescale_clean=True,
                           do_revert_to_singlescale=True,
                           do_singlescale_mask=True,
                           singlescale_mask_absolute=absolute,
                           singlescale_mask_high_snr=high_snr,
                           singlescale_mask_low_snr=low_snr,
                           convergence_fracflux=convergence_fracflux,
                           singlescale_threshold_value=singlescale_threshold_value,
                           do_recombine_cubes=True,
                           do_export_to_fits=True,
                           do_cleanup=True,
                           overwrite=True,
                           )

im_handler = ImagingChunkedHandler(
    target=target,
    config=interf_config,
    product=line_product,
    key_handler=key_handler,
    imaging_method=imaging_method,
    chunksize=chunksize,
    make_temp_dir=True,
    temp_path=chunk_temp_path,
    copy_ms_to_temp=False,
)

im_handler.task_complete_gather_into_cubes()
im_handler.task_export_to_fits()
im_handler.task_cleanup()

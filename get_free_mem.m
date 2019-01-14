function mem=get_free_mem()
    [~,out]=system('vm_stat | grep "Pages free"');
    mem=sscanf(out,'Pages free: %f.');
    mem=mem*4096/1000000000;
end
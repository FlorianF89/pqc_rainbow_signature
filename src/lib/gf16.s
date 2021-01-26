	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 10, 15	sdk_version 10, 15, 6
	.globl	_bitsliced_gf16_is_one  ## -- Begin function bitsliced_gf16_is_one
	.p2align	4, 0x90
_bitsliced_gf16_is_one:                 ## @bitsliced_gf16_is_one
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	cmpq	$1, 16(%rbp)
	jne	LBB0_5
## %bb.1:
	leaq	16(%rbp), %rax
	cmpq	$0, 24(%rax)
	jne	LBB0_5
## %bb.2:
	cmpq	$0, 8(%rax)
	jne	LBB0_5
## %bb.3:
	cmpq	$0, 16(%rax)
	je	LBB0_4
LBB0_5:
	xorl	%eax, %eax
	popq	%rbp
	retq
LBB0_4:
	movl	$1, %eax
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_gf16_is_zero           ## -- Begin function gf16_is_zero
	.p2align	4, 0x90
_gf16_is_zero:                          ## @gf16_is_zero
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	24(%rbp), %rcx
	orq	16(%rbp), %rcx
	orq	32(%rbp), %rcx
	orq	40(%rbp), %rcx
                                        ## kill: def $edi killed $edi def $rdi
	xorl	%eax, %eax
	btq	%rdi, %rcx
	setae	%al
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_copy_gf16              ## -- Begin function copy_gf16
	.p2align	4, 0x90
_copy_gf16:                             ## @copy_gf16
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	(%rsi), %rax
	movq	%rax, (%rdi)
	movq	8(%rsi), %rax
	movq	%rax, 8(%rdi)
	movq	16(%rsi), %rax
	movq	%rax, 16(%rdi)
	movq	24(%rsi), %rax
	movq	%rax, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_move_two_halves_gf16_into_one ## -- Begin function move_two_halves_gf16_into_one
	.p2align	4, 0x90
_move_two_halves_gf16_into_one:         ## @move_two_halves_gf16_into_one
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	(%rdx), %rax
	shlq	$32, %rax
	orq	(%rsi), %rax
	movq	%rax, (%rdi)
	movq	8(%rdx), %rax
	shlq	$32, %rax
	orq	8(%rsi), %rax
	movq	%rax, 8(%rdi)
	movq	16(%rdx), %rax
	shlq	$32, %rax
	orq	16(%rsi), %rax
	movq	%rax, 16(%rdi)
	movq	24(%rdx), %rax
	shlq	$32, %rax
	orq	24(%rsi), %rax
	movq	%rax, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_bitsliced_addition     ## -- Begin function bitsliced_addition
	.p2align	4, 0x90
_bitsliced_addition:                    ## @bitsliced_addition
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	(%rdx), %rax
	xorq	(%rsi), %rax
	movq	%rax, (%rdi)
	movq	8(%rdx), %rax
	xorq	8(%rsi), %rax
	movq	%rax, 8(%rdi)
	movq	16(%rdx), %rax
	xorq	16(%rsi), %rax
	movq	%rax, 16(%rdi)
	movq	24(%rdx), %rax
	xorq	24(%rsi), %rax
	movq	%rax, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_shift_left_gf16        ## -- Begin function shift_left_gf16
	.p2align	4, 0x90
_shift_left_gf16:                       ## @shift_left_gf16
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movl	%edx, %ecx
	movq	(%rsi), %rax
	shlq	%cl, %rax
	movq	%rax, (%rdi)
	movq	8(%rsi), %rax
	shlq	%cl, %rax
	movq	%rax, 8(%rdi)
	movq	16(%rsi), %rax
	shlq	%cl, %rax
	movq	%rax, 16(%rdi)
	movq	24(%rsi), %rax
                                        ## kill: def $cl killed $cl killed $ecx
	shlq	%cl, %rax
	movq	%rax, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_shift_right_gf16       ## -- Begin function shift_right_gf16
	.p2align	4, 0x90
_shift_right_gf16:                      ## @shift_right_gf16
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movl	%edx, %ecx
	movq	(%rsi), %rax
	shrq	%cl, %rax
	movq	%rax, (%rdi)
	movq	8(%rsi), %rax
	shrq	%cl, %rax
	movq	%rax, 8(%rdi)
	movq	16(%rsi), %rax
	shrq	%cl, %rax
	movq	%rax, 16(%rdi)
	movq	24(%rsi), %rax
                                        ## kill: def $cl killed $cl killed $ecx
	shrq	%cl, %rax
	movq	%rax, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_extract_one_gf16_element_and_place_it_in_given_position ## -- Begin function extract_one_gf16_element_and_place_it_in_given_position
	.p2align	4, 0x90
_extract_one_gf16_element_and_place_it_in_given_position: ## @extract_one_gf16_element_and_place_it_in_given_position
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movl	%ecx, %eax
	movl	%esi, %r8d
	movq	(%rdx), %rsi
	shrq	%cl, %rsi
	andl	$1, %esi
	movl	%r8d, %ecx
	shlq	%cl, %rsi
	movq	%rsi, (%rdi)
	movq	8(%rdx), %rsi
	movl	%eax, %ecx
	shrq	%cl, %rsi
	andl	$1, %esi
	movl	%r8d, %ecx
	shlq	%cl, %rsi
	movq	%rsi, 8(%rdi)
	movq	16(%rdx), %rsi
	movl	%eax, %ecx
	shrq	%cl, %rsi
	andl	$1, %esi
	movl	%r8d, %ecx
	shlq	%cl, %rsi
	movq	%rsi, 16(%rdi)
	movq	24(%rdx), %rdx
	movl	%eax, %ecx
	shrq	%cl, %rdx
	andl	$1, %edx
	movl	%r8d, %ecx
	shlq	%cl, %rdx
	movq	%rdx, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_bitsliced_multiplication ## -- Begin function bitsliced_multiplication
	.p2align	4, 0x90
_bitsliced_multiplication:              ## @bitsliced_multiplication
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r14
	pushq	%rbx
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	movq	(%rdx), %r8
	andq	(%rsi), %r8
	movq	%r8, (%rdi)
	movq	(%rsi), %r10
	movq	8(%rdx), %rcx
	andq	%r10, %rcx
	movq	(%rdx), %rax
	movq	8(%rsi), %r9
	andq	%rax, %r9
	xorq	%rcx, %r9
	movq	%r9, 8(%rdi)
	movq	16(%rdx), %rcx
	andq	%r10, %rcx
	movq	16(%rsi), %r11
	andq	%rax, %r11
	xorq	%rcx, %r11
	movq	%r11, 16(%rdi)
	andq	24(%rdx), %r10
	movq	8(%rsi), %r14
	movq	16(%rdx), %rbx
	andq	%r14, %rbx
	xorq	%r10, %rbx
	movq	8(%rdx), %r10
	movq	16(%rsi), %rcx
	andq	%r10, %rcx
	andq	24(%rsi), %rax
	xorq	%rcx, %rax
	xorq	%rbx, %rax
	movq	%rax, 24(%rdi)
	andq	%r14, %r10
	xorq	%r10, %r8
	xorq	%r9, %r10
	movq	%r10, 8(%rdi)
	movq	24(%rdx), %rbx
	andq	8(%rsi), %rbx
	xorq	%rbx, %r11
	movq	%r11, 16(%rdi)
	xorq	%rax, %rbx
	movq	%rbx, 24(%rdi)
	movq	16(%rdx), %rcx
	andq	16(%rsi), %rcx
	xorq	%rcx, %r11
	movq	%r11, 16(%rdi)
	movq	24(%rdx), %rax
	andq	16(%rsi), %rax
	xorq	%rax, %rcx
	xorq	%r10, %rcx
	movq	%rcx, 8(%rdi)
	xorq	%rax, %rbx
	movq	%rbx, 24(%rdi)
	movq	8(%rdx), %r9
	andq	24(%rsi), %r9
	xorq	%r9, %r11
	movq	%r11, 16(%rdi)
	xorq	%rbx, %r9
	movq	%r9, 24(%rdi)
	movq	16(%rdx), %rbx
	andq	24(%rsi), %rbx
	xorq	%rbx, %rax
	xorq	%r8, %rax
	xorq	%rbx, %rcx
	movq	%rcx, 8(%rdi)
	xorq	%r9, %rbx
	movq	%rbx, 24(%rdi)
	movq	24(%rdx), %rcx
	andq	24(%rsi), %rcx
	xorq	%rcx, %rax
	movq	%rax, (%rdi)
	xorq	%rcx, %r11
	movq	%r11, 16(%rdi)
	xorq	%rbx, %rcx
	movq	%rcx, 24(%rdi)
	popq	%rbx
	popq	%r14
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_bitsliced_muladd       ## -- Begin function bitsliced_muladd
	.p2align	4, 0x90
_bitsliced_muladd:                      ## @bitsliced_muladd
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%rbx
	.cfi_offset %rbx, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movq	(%rdx), %r8
	andq	(%rsi), %r8
	xorq	(%rdi), %r8
	movq	%r8, (%rdi)
	movq	(%rsi), %r9
	movq	8(%rdx), %rcx
	andq	%r9, %rcx
	movq	(%rdx), %rax
	movq	8(%rsi), %r10
	andq	%rax, %r10
	xorq	%rcx, %r10
	xorq	8(%rdi), %r10
	movq	%r10, 8(%rdi)
	movq	16(%rdx), %r11
	andq	%r9, %r11
	movq	16(%rsi), %r15
	andq	%rax, %r15
	xorq	%r11, %r15
	xorq	16(%rdi), %r15
	movq	%r15, 16(%rdi)
	andq	24(%rdx), %r9
	movq	8(%rsi), %r11
	movq	16(%rdx), %rbx
	andq	%r11, %rbx
	xorq	%r9, %rbx
	movq	8(%rdx), %r9
	movq	16(%rsi), %r14
	andq	%r9, %r14
	andq	24(%rsi), %rax
	xorq	%r14, %rax
	xorq	%rbx, %rax
	xorq	24(%rdi), %rax
	movq	%rax, 24(%rdi)
	andq	%r11, %r9
	xorq	%r9, %r10
	movq	%r10, 8(%rdi)
	movq	24(%rdx), %r11
	andq	8(%rsi), %r11
	xorq	%r11, %r15
	movq	%r15, 16(%rdi)
	xorq	%rax, %r11
	movq	%r11, 24(%rdi)
	movq	16(%rdx), %rax
	andq	16(%rsi), %rax
	xorq	%rax, %r15
	movq	%r15, 16(%rdi)
	movq	24(%rdx), %rbx
	andq	16(%rsi), %rbx
	xorq	%rbx, %r9
	xorq	%rbx, %rax
	xorq	%r10, %rax
	movq	%rax, 8(%rdi)
	xorq	%r11, %rbx
	movq	%rbx, 24(%rdi)
	movq	8(%rdx), %rcx
	andq	24(%rsi), %rcx
	xorq	%rcx, %r15
	movq	%r15, 16(%rdi)
	xorq	%rbx, %rcx
	movq	%rcx, 24(%rdi)
	movq	16(%rdx), %rbx
	andq	24(%rsi), %rbx
	xorq	%rbx, %rax
	movq	%rax, 8(%rdi)
	xorq	%rbx, %rcx
	movq	%rcx, 24(%rdi)
	movq	24(%rdx), %rax
	andq	24(%rsi), %rax
	xorq	%r8, %r9
	xorq	%rax, %rbx
	xorq	%r9, %rbx
	movq	%rbx, (%rdi)
	xorq	%rax, %r15
	movq	%r15, 16(%rdi)
	xorq	%rcx, %rax
	movq	%rax, 24(%rdi)
	popq	%rbx
	popq	%r14
	popq	%r15
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	__gf16v_madd_u32        ## -- Begin function _gf16v_madd_u32
	.p2align	4, 0x90
__gf16v_madd_u32:                       ## @_gf16v_madd_u32
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$56, %rsp
	.cfi_offset %rbx, -56
	.cfi_offset %r12, -48
	.cfi_offset %r13, -40
	.cfi_offset %r14, -32
	.cfi_offset %r15, -24
	movl	%ecx, %r15d
	movl	%edx, %r12d
	movq	%rsi, %r14
	movq	%rdi, %r13
	leaq	L_.str(%rip), %rdi
	movl	%edx, %esi
	movl	%ecx, %edx
	xorl	%eax, %eax
	callq	_printf
	## InlineAsm Start
	.byte	15
	.byte	49
	## InlineAsm End
	movq	%rax, -48(%rbp)         ## 8-byte Spill
	addl	$31, %r15d
	shrl	$5, %r15d
	je	LBB10_3
## %bb.1:
	movl	%r12d, %eax
	andl	$1, %eax
	negl	%eax
	cltq
	movq	%rax, -88(%rbp)         ## 8-byte Spill
	movl	%r12d, %eax
	shrl	%eax
	andl	$1, %eax
	negl	%eax
	cltq
	movq	%rax, -80(%rbp)         ## 8-byte Spill
	movl	%r12d, %eax
	shrl	$2, %eax
	andl	$1, %eax
	negl	%eax
	cltq
	movq	%rax, -72(%rbp)         ## 8-byte Spill
	shrl	$3, %r12d
	andl	$1, %r12d
	negl	%r12d
	movslq	%r12d, %rax
	movq	%rax, -64(%rbp)         ## 8-byte Spill
	movl	%r15d, %eax
	shlq	$5, %rax
	movq	%rax, -56(%rbp)         ## 8-byte Spill
	xorl	%ecx, %ecx
	movq	-80(%rbp), %r11         ## 8-byte Reload
	movq	-64(%rbp), %r12         ## 8-byte Reload
	.p2align	4, 0x90
LBB10_2:                                ## =>This Inner Loop Header: Depth=1
	movq	(%r14,%rcx), %r15
	movq	-88(%rbp), %rbx         ## 8-byte Reload
	andq	%rbx, %r15
	xorq	(%r13,%rcx), %r15
	movq	%r15, (%r13,%rcx)
	movq	(%r14,%rcx), %rdi
	movq	%rdi, %rax
	movq	-72(%rbp), %rdx         ## 8-byte Reload
	andq	%rdx, %rax
	movq	%rdx, %r10
	movq	8(%r14,%rcx), %rdx
	andq	%rbx, %rdx
	xorq	%rax, %rdx
	xorq	8(%r13,%rcx), %rdx
	movq	%rdx, 8(%r13,%rcx)
	movq	%rdi, %rsi
	andq	%r11, %rsi
	movq	16(%r14,%rcx), %rax
	andq	%rbx, %rax
	xorq	%rsi, %rax
	xorq	16(%r13,%rcx), %rax
	movq	%rax, 16(%r13,%rcx)
	andq	%r12, %rdi
	movq	8(%r14,%rcx), %r9
	movq	%r9, %rsi
	andq	%r11, %rsi
	xorq	%rdi, %rsi
	movq	16(%r14,%rcx), %rdi
	andq	%r10, %rdi
	movq	24(%r14,%rcx), %r8
	andq	%rbx, %r8
	xorq	%rdi, %r8
	xorq	%rsi, %r8
	xorq	24(%r13,%rcx), %r8
	andq	%r10, %r9
	movq	%r10, %rbx
	xorq	%r9, %rdx
	movq	%rdx, 8(%r13,%rcx)
	movq	8(%r14,%rcx), %rsi
	andq	%r12, %rsi
	xorq	%rsi, %rax
	movq	%rax, 16(%r13,%rcx)
	movq	16(%r14,%rcx), %r10
	andq	%r11, %r10
	xorq	%r10, %rax
	movq	%rax, 16(%r13,%rcx)
	movq	16(%r14,%rcx), %rdi
	andq	%r12, %rdi
	xorq	%rdi, %r9
	xorq	%r15, %r9
	xorq	%rdi, %r10
	xorq	%rsi, %rdi
	xorq	%r8, %rdi
	movq	%rdi, 24(%r13,%rcx)
	movq	24(%r14,%rcx), %rsi
	andq	%rbx, %rsi
	xorq	%rsi, %rdi
	movq	%rdi, 24(%r13,%rcx)
	movq	24(%r14,%rcx), %rbx
	andq	%r11, %rbx
	xorq	%rbx, %r10
	xorq	%rdx, %r10
	movq	%r10, 8(%r13,%rcx)
	xorq	%rbx, %rdi
	movq	%rdi, 24(%r13,%rcx)
	movq	24(%r14,%rcx), %rdx
	andq	%r12, %rdx
	xorq	%rdx, %rbx
	xorq	%r9, %rbx
	movq	%rbx, (%r13,%rcx)
	xorq	%rdx, %rsi
	xorq	%rax, %rsi
	movq	%rsi, 16(%r13,%rcx)
	xorq	%rdi, %rdx
	movq	%rdx, 24(%r13,%rcx)
	addq	$32, %rcx
	cmpq	%rcx, -56(%rbp)         ## 8-byte Folded Reload
	jne	LBB10_2
LBB10_3:
	## InlineAsm Start
	.byte	15
	.byte	49
	## InlineAsm End
	subq	-48(%rbp), %rax         ## 8-byte Folded Reload
	leaq	L_.str.1(%rip), %rdi
	movq	%rax, %rsi
	xorl	%eax, %eax
	addq	$56, %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	jmp	_printf                 ## TAILCALL
	.cfi_endproc
                                        ## -- End function
	.globl	_bitsliced_vectorized_multiplication ## -- Begin function bitsliced_vectorized_multiplication
	.p2align	4, 0x90
_bitsliced_vectorized_multiplication:   ## @bitsliced_vectorized_multiplication
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_bitsliced_square       ## -- Begin function bitsliced_square
	.p2align	4, 0x90
_bitsliced_square:                      ## @bitsliced_square
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	movq	8(%rsi), %rax
	movq	24(%rsi), %rcx
	movq	(%rsi), %rdx
	xorq	%rcx, %rdx
	xorq	%rax, %rdx
	movq	%rdx, (%rdi)
	movq	16(%rsi), %rdx
	xorq	%rdx, %rax
	movq	%rax, 8(%rdi)
	xorq	%rcx, %rdx
	movq	%rdx, 16(%rdi)
	movq	%rcx, 24(%rdi)
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_bitsliced_inversion    ## -- Begin function bitsliced_inversion
	.p2align	4, 0x90
_bitsliced_inversion:                   ## @bitsliced_inversion
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r14
	pushq	%rbx
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	movq	16(%rsi), %r10
	movq	24(%rsi), %r8
	movq	8(%rsi), %rax
	xorq	%r10, %rax
	movq	(%rsi), %rcx
	xorq	%rax, %rcx
	xorq	%r8, %rax
	movq	%rax, 8(%rdi)
	movq	%r10, 16(%rdi)
	movq	16(%rsi), %r9
	xorq	%r9, %r8
	movq	%r8, 24(%rdi)
	movq	24(%rsi), %rdx
	movq	8(%rsi), %r11
	andq	%rdx, %r11
	xorq	%rcx, %r11
	movq	%r11, (%rdi)
	movq	(%rsi), %rcx
	andq	%rcx, %r9
	xorq	%rax, %r9
	movq	%r9, 8(%rdi)
	andq	%rdx, %rcx
	xorq	%rcx, %r10
	movq	%r10, 16(%rdi)
	xorq	%rcx, %r8
	movq	%r8, 24(%rdi)
	movq	8(%rsi), %rbx
	movq	16(%rsi), %r14
	andq	%rbx, %r14
	movq	%rcx, %rax
	xorq	%r14, %rax
	xorq	%r11, %rax
	xorq	%r14, %r10
	movq	%r10, 16(%rdi)
	andq	%rcx, %rbx
	andq	16(%rsi), %rcx
	xorq	%rbx, %rax
	movq	24(%rsi), %rdx
	andq	%r14, %rdx
	xorq	%rdx, %rax
	xorq	%rcx, %rax
	movq	%rax, (%rdi)
	xorq	%rdx, %rbx
	xorq	%r9, %rbx
	movq	%rbx, 8(%rdi)
	xorq	%rdx, %r10
	xorq	%rcx, %r10
	movq	%r10, 16(%rdi)
	xorq	%r8, %rdx
	movq	%rdx, 24(%rdi)
	andq	(%rsi), %r14
	xorq	%rax, %r14
	movq	%r14, (%rdi)
	popq	%rbx
	popq	%r14
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.globl	_main                   ## -- Begin function main
	.p2align	4, 0x90
_main:                                  ## @main
	.cfi_startproc
## %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset %rbp, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register %rbp
	pushq	%r14
	pushq	%rbx
	subq	$64, %rsp
	.cfi_offset %rbx, -32
	.cfi_offset %r14, -24
	leaq	L_.str.2(%rip), %rdi
	leaq	L_.str.3(%rip), %rsi
	callq	_fopen
	movq	%rax, %rbx
	leaq	-48(%rbp), %r14
	movl	$1, %esi
	movl	$8, %edx
	movq	%r14, %rdi
	movq	%rax, %rcx
	callq	_fread
	leaq	-32(%rbp), %rdi
	movl	$1, %esi
	movl	$8, %edx
	movq	%rbx, %rcx
	callq	_fread
	leaq	-40(%rbp), %rdi
	movl	$1, %esi
	movl	$8, %edx
	movq	%rbx, %rcx
	callq	_fread
	leaq	-24(%rbp), %rdi
	movl	$1, %esi
	movl	$8, %edx
	movq	%rbx, %rcx
	callq	_fread
	movq	%rbx, %rdi
	callq	_fclose
	leaq	-80(%rbp), %rbx
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$16, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	movq	%rbx, %rdi
	movq	%r14, %rsi
	xorl	%edx, %edx
	movl	$32, %ecx
	callq	__gf16v_madd_u32
	xorl	%eax, %eax
	addq	$64, %rsp
	popq	%rbx
	popq	%r14
	popq	%rbp
	retq
	.cfi_endproc
                                        ## -- End function
	.section	__TEXT,__cstring,cstring_literals
L_.str:                                 ## @.str
	.asciz	"_gf16v_madd_u32 got used with %02X multiplying a %d byte vector: "

L_.str.1:                               ## @.str.1
	.asciz	"Time is %d\n"

L_.str.2:                               ## @.str.2
	.asciz	"/dev/urandom"

L_.str.3:                               ## @.str.3
	.asciz	"rb"

.subsections_via_symbols
